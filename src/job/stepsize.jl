"""
    stepsize(configfile)

Experiments on the stepsize extension of DSEA, with constant, decaying, and adaptive step
size strategies.
"""
function stepsize(configfile::String)
    
    # read configuration
    c = load_file(configfile)
    B = c["num_bootstraps"]
    K = c["dsea_iterations"]
    outfile     = c["outfile"]
    experiments = c["experiments"]
    
    # generate seeds for experiments and bootstraps
    srand(c["seed"])
    pushseeds!(experiments, B)
    for (i_exp, exp) in enumerate(experiments) # also add an info string
        exp["info"] = "$(i_exp)/$(length(experiments))"
    end
    
    # read data
    dataset = Data.dataset(c["dataset"]["id"], nobs = -1, nobs_train = -1) # full data set
    discr   = Data.discretizer(dataset)
    bins    = Data.bins(discr)
    X_full  = Data.X_data(dataset)
    y_full  = encode(discr, Data.y_data(dataset))
    
    X_auxtrain, y_auxtrain = try
        (Data.X_train(dataset), encode(discr, Data.y_train(dataset)))
    catch (nothing, nothing) end
    
    # more data-dependent configuration
    skconfig     = joinpath("conf/skl", c["dataset"]["skconfig"] * ".yml")
    f_train_conf = c["dataset"]["f_train"]
    
    # split into training and observed data sets
    X_data, y_data, X_train, y_train = Util.shuffle_split_subsample( X_full, y_full,
                                                                     X_auxtrain, y_auxtrain,
                                                                     f_train_conf )
    
    # parallel execution
    info("Starting $(length(experiments)) experiments on $(nworkers()) worker(s).")
    df = vcat(pmap(exp -> _stepsize_job(exp, X_data, y_data, X_train, y_train, bins,
                                        f_train_conf, K, skconfig),
                   experiments)...)
                 # on_error = err -> isa(err, InterruptException) ? rethrow(err) : DataFrame(err = string(err), trace = string(catch_stacktrace()))
    
    # output
    writetable(outfile, df)
    info("Results written to $outfile")
    df
    
end

# independent sub-routine of stepsize
function _stepsize_job{TN<:Number,TI<:Int}( exp      :: Dict{Any,Any},
                                            X_data   :: Matrix{TN},
                                            y_data   :: Array{TI,1}, 
                                            X_train  :: Matrix{TN},
                                            y_train  :: Array{TI,1},
                                            bins     :: AbstractArray{TI,1},
                                            f_train  :: String,
                                            K        :: Int,
                                            skconfig :: String )
    
    # read basic configuration
    srand(exp["seed"])
    seeds           = exp["bootstrap_seeds"]
    fixweighting    = exp["fixweighting"]
    stepsize_conf   = exp["stepsize"] # nested config
    stepsize_method = stepsize_conf["method"]
    
    # informative name of strategy requires configuration
    name = exp["name"]
    if stepsize_method == "constant"
        name = replace(name, "\$alpha", stepsize_conf["alpha"])
    elseif stepsize_method == "decay_mul" || stepsize_method == "decay_exp"
        name = replace(name, "\$eta", stepsize_conf["eta"])
    elseif stepsize_method == "run"
        name = replace(name, "\$tau", Util.latex_e(stepsize_conf["tau"], dollars = false))
    end
    
    # bootstrap
    progress = DataFrame( name         = String[],
                          method       = String[],
                          fixweighting = Bool[],
                          b            = Int64[],
                          k            = Int64[],
                          alpha        = Float64[],
                          f            = Array{Float64, 1}[] )
    for b in 1:length(seeds)
        info("Starting bootstrap $b/$(length(seeds)) of experiment $(exp["info"]) on $name")
        rng = MersenneTwister(seeds[b])
        itrain, _   = Util.bootstrap_sample_indices(rng, length(y_train))
        idata,  oob = Util.bootstrap_sample_indices(rng, length(y_data))
        
        # true spectrum in this bootstrap sample (in-bag and oob)
        f_ib = fit_pdf(y_data[idata], bins) # true in-bag spectrum is re-used for optimal step
        for (truth_name, indices) in zip([TRUESPEC_IB, TRUESPEC_OOB], [idata, oob])
            f_true = truth_name == TRUESPEC_IB ? f_ib : fit_pdf(y_data[indices], bins)
            push!(progress, [ truth_name, "", false, b, -1, NaN, f_true ])
        end
        
        # line search strategy
        alpha = _stepsize_alpha(stepsize_conf, X_data[idata, :], X_train[itrain, :],
                                y_train[itrain], bins, f_ib)
        
        # configure classifier and set up inspection
        tp = train_and_predict_proba(Util.classifier_from_config(skconfig)) # from CherenkovDeconvolution.Sklearn
        inspect = (f, k, chi2s, alpha) -> push!(progress, [ name, stepsize_method, fixweighting, b, k, alpha, f ])
        
        # 
        # DSEA
        # 
        # Deconvolve the configured data set to investigate the benefit of the line search
        # extension to DSEA.
        # 
        dsea(X_data[idata, :], X_train[itrain, :], y_train[itrain], tp, bins,
             alpha = alpha, fixweighting = fixweighting, K = K, inspect = inspect)
    end
    return progress
    
end

# return the alpha argument for DSEA, which may be a constant value or a function
_stepsize_alpha{TN<:Number,TI<:Int}( stepsize_conf :: Dict{Any, Any},
                                     X_data        :: Matrix{TN},
                                     X_train       :: Matrix{TN},
                                     y_train       :: Array{TI,1},
                                     bins          :: AbstractArray{TI,1},
                                     f_ib          :: Array{Float64,1} ) =
    if stepsize_conf["method"] == "constant"
        stepsize_conf["alpha"] # fixed value
        
    elseif stepsize_conf["method"] == "decay_mul"
        alpha_decay_mul(stepsize_conf["eta"], stepsize_conf["a_1"]) # slow decay function
        
    elseif stepsize_conf["method"] == "decay_exp"
        alpha_decay_exp(stepsize_conf["eta"], stepsize_conf["a_1"]) # fast decay function
        
    elseif stepsize_conf["method"] == "optimal" # minimizer of true(!) in-bag error
        (k::Int, pk::Array{Float64,1}, f::Array{Float64,1}) -> begin
            a_min, a_max = CherenkovDeconvolution._alpha_range(pk, f)
            Optim.optimize(a -> _stepsize_emd(f_ib, f + a * pk), a_min, a_max).minimizer
        end
    
    elseif stepsize_conf["method"] == "run" # maximize likelihood (like RUN) in DSEA direction
        
        # discretize the feature space
        discr_x = TreeDiscretizer(X_train, y_train, stepsize_conf["num_clusters"])
        bins_x  = Sklearn.bins(discr_x)
        x_data  = encode(discr_x, X_data)
        x_train = encode(discr_x, X_train)
        
        # return a function object from CherenkovDeconvolution.jl
        alpha_adaptive_run(x_data, x_train, y_train, bins, stepsize_conf["tau"])
        
    else
        error("Stepsize method '$method' is not available")
    end

# emd between two un-normalized estimates
_stepsize_emd(a::Array{Float64,1}, b::Array{Float64,1}) = Util.mdpa(normalizepdf(a, b)...)


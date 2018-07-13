"""
    smoothing(configfile="conf/job/smoothing.yml")

Experiments on the smoothing extension of DSEA.
"""
function smoothing(configfile::String="conf/job/smoothing.yml")
    
    # read configuration
    c = load_file(configfile)
    B = c["num_bootstraps"]
    K = c["dsea_iterations"]
    outfile = c["outfile"]
    
    # expand experiments
    experiments = vcat(map(
        exp -> begin
            if exp["smoothing"]["method"] == "polynomial"
                Yaml.expand(exp, ["smoothing", "order"])
            else # e.g., no smoothing
                exp
            end
        end, c["experiments"])...)
    
    # generate seeds for experiments and bootstraps
    srand(c["seed"])
    pushseeds!(experiments, B)
    for (i_exp, exp) in enumerate(experiments) # also add an info string
        exp["info"] = "$(i_exp)/$(length(experiments))"
    end
    
    # read data
    dataset = Data.Gaussian()
    discr   = Data.discretizer(dataset)
    bins    = Data.bins(discr)
    X_data  = Data.X_data(dataset)
    X_train = Data.X_train(dataset)
    y_data  = encode(discr, Data.y_data(dataset))
    y_train = encode(discr, Data.y_train(dataset))
    
    # parallel execution
    info("Starting $(length(experiments)) experiments on $(nworkers()) worker(s).")
    df = vcat(pmap(exp -> _smoothing_job(exp, X_data, y_data, X_train, y_train, bins, K),
                   experiments)...)
    
    # output
    writetable(outfile, df)
    info("Results written to $outfile")
    df
    
end

# independent sub-routine of smoothing
function _smoothing_job{TN<:Number,TI<:Int}( exp      :: Dict{Any,Any},
                                             X_data   :: Matrix{TN},
                                             y_data   :: Array{TI,1}, 
                                             X_train  :: Matrix{TN},
                                             y_train  :: Array{TI,1},
                                             bins     :: AbstractArray{TI,1},
                                             K        :: Int )
    # read basic configuration
    srand(exp["seed"])
    seeds            = exp["bootstrap_seeds"]
    fixweighting     = exp["fixweighting"]
    smoothing_conf   = exp["smoothing"]
    smoothing_method = smoothing_conf["method"]
    
    # informative name of strategy requires configuration
    name = exp["name"]
    if smoothing_method == "polynomial"
        name = replace(name, "\$order", smoothing_conf["order"])
    end
    
    # prepare smoothing function
    smoothing_func = if smoothing_method == "polynomial"
        polynomial_smoothing(smoothing_conf["order"]) # from CherenkovDeconvolution.Util
    elseif smoothing_method == "none"
        Base.identity
    else
        throw(ArgumentError("Smoothing method '$(smoothing_conf["method"])' not available"))
    end
    
    # bootstrap
    progress = DataFrame( name         = String[],
                          method       = String[],
                          fixweighting = Bool[],
                          b            = Int64[],
                          k            = Int64[],
                          f            = Array{Float64, 1}[] )
    for b in 1:length(seeds)
        info("Starting bootstrap $b/$(length(seeds)) of experiment $(exp["info"]) on '$name' smoothing")
        rng = MersenneTwister(seeds[b])
        itrain, _   = Util.bootstrap_sample_indices(rng, length(y_train))
        idata,  oob = Util.bootstrap_sample_indices(rng, length(y_data))
        
        # true spectrum in this bootstrap sample (in-bag and oob)
        for (truth_name, indices) in zip([TRUESPEC_IB, TRUESPEC_OOB], [idata, oob])
            f_true = fit_pdf(y_data[indices], bins)
            push!(progress, [ truth_name, "", false, b, -1, f_true ])
        end
        
        # configure classifier and set up inspection
        tp = train_and_predict_proba(Util.classifier_from_config("conf/skl/nb.yml"))
        inspect = (f, k, chi2s, alpha) -> push!(progress, [ name, smoothing_method, fixweighting, b, k, f ])
                                                            
        # deconvolve the configured data set, representing one amount of smearing.
        dsea(X_data[idata, :], X_train[itrain, :], y_train[itrain], tp, bins,
             smoothing = smoothing_func, fixweighting = fixweighting, K = K, inspect = inspect)
    end
    return progress
    
end


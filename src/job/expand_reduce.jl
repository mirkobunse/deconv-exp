"""
    expand_reduce(configfile="conf/job/expand_reduce.yml")

Experiments on expansion/reduction in DSEA.
"""
function expand_reduce(configfile::String="conf/job/expand_reduce.yml")
    
    # read configuration
    c = load_file(configfile)
    B = c["num_bootstraps"]
    K = c["dsea_iterations"]
    experiments = c["experiments"]
    outfile = c["outfile"]
    
    # generate seeds for experiments and bootstraps
    srand(c["seed"])
    pushseeds!(experiments, B)
    for (i_exp, exp) in enumerate(experiments) # also add an info string
        exp["info"] = "$(i_exp)/$(length(experiments))"
    end
    
    # read data
    dataset   = Data.Gaussian()
    discr     = Data.discretizer(dataset)
    X_data    = Data.X_data(dataset)
    X_train   = Data.X_train(dataset)
    y_data_c  = Data.y_data(dataset) # continuous values
    y_train_c = Data.y_train(dataset)
    
    # parallel execution
    info("Starting $(length(experiments)) experiments on $(nworkers()) worker(s).")
    df = vcat(pmap(exp -> _expand_reduce_job(exp, X_data, y_data_c, X_train, y_train_c, discr, K),
                   experiments)...)
    
    # output
    writetable(outfile, df)
    info("Results written to $outfile")
    df
    
end

# independent sub-routine of smoothing
function _expand_reduce_job{TN<:Number,TL<:Number}( exp       :: Dict{Any,Any},
                                                    X_data    :: Matrix{TN},
                                                    y_data_c  :: Array{TL,1}, 
                                                    X_train   :: Matrix{TN},
                                                    y_train_c :: Array{TL,1},
                                                    discr     :: AbstractDiscretizer,
                                                    K         :: Int )
    # read basic configuration
    srand(exp["seed"])
    seeds  = exp["bootstrap_seeds"]
    factor = exp["factor"]
    name   = exp["name"]
    
    # discretize truth as usual
    y_data = encode(discr, y_data_c)
    bins_f = Data.bins(discr) # bins in original solution space
    
    # discretize y_train, considering the expand/reduce configuration
    discr   = expansion_discretizer(discr, factor)
    y_train = encode(discr, y_train_c)
    bins    = Data.bins(discr) # bins in expansion
    
    # bootstrap
    progress = DataFrame( name   = String[],
                          factor = Int64[],
                          b      = Int64[],
                          k      = Int64[],
                          f      = Array{Float64, 1}[] )
    for b in 1:length(seeds)
        info("Starting bootstrap $b/$(length(seeds)) of experiment $(exp["info"]) on '$name'")
        rng = MersenneTwister(seeds[b])
        itrain, _   = Util.bootstrap_sample_indices(rng, length(y_train))
        idata,  oob = Util.bootstrap_sample_indices(rng, length(y_data))
        
        # true spectrum in this bootstrap sample (in-bag and oob)
        for (truth_name, indices) in zip([TRUESPEC_IB, TRUESPEC_OOB], [idata, oob])
            f_true = fit_pdf(y_data[indices], bins_f)
            push!(progress, [ truth_name, -1, b, -1, f_true ])
        end
        
        # line search strategy
        stepsize_conf = Dict{Any,Any}( "method"       => "run",
                                       "num_clusters" => length(bins),
                                       "tau"          => 0.0 )
        alpha = _stepsize_alpha(stepsize_conf, X_data[idata, :], X_train[itrain, :],
                                y_train[itrain], bins, zeros(0))
        
        # configure classifier and wrap usual inspection for expansion
        tp = train_and_predict_proba(Util.classifier_from_config("conf/skl/nb.yml"))
        inspect_usual = (f, k, chi2s, alpha) -> push!(progress, [ name, factor, b, k, f ])
        inspect       = inspect_expansion(inspect_usual, factor) # from CherenkovDeconvolution.jl
                                                            
        # deconvolve the configured data set, representing one amount of smearing.
        dsea(X_data[idata, :], X_train[itrain, :], y_train[itrain], tp, bins,
             alpha = alpha, K = K, inspect = inspect)
    end
    return progress
    
end


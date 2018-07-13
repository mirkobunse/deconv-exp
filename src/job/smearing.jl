"""
    smearing(configfile="conf/job/smearing.yml")

Basic DSEA experiment for data with different amounts of smearing.
"""
function smearing(configfile::String="conf/job/smearing.yml")
    
    # configuration
    c = load_file(configfile)
    B = c["num_bootstraps"]
    K = c["dsea_iterations"]
    outfile     = c["outfile"]
    experiments = c["experiments"]
    
    # generate seeds for experiments and bootstrap
    srand(c["seed"])
    pushseeds!(experiments, B)
    
    # parallel execution
    info("Running $(length(experiments)) experiments independently on $(nworkers()) worker(s)")
    df = vcat(pmap(exp -> _smearing_job(exp, K), experiments)...)
    
    # output
    writetable(outfile, df)
    info("Results written to $outfile")
    df
    
end

# independent sub-routine of smearing
function _smearing_job(exp::Dict{Any,Any}, K::Int)
    srand(exp["seed"])
    
    # read data, discretizing the target variable
    configfile = exp["configfile"]
    name       = exp["name"]
    dataset = Data.Gaussian(configfile)
    discr   = Data.discretizer(dataset)
    bins    = Data.bins(discr)
    X_data  = Data.X_data(dataset)
    X_train = Data.X_train(dataset)
    y_data  = encode(discr, Data.y_data(dataset))
    y_train = encode(discr, Data.y_train(dataset))
    
    # bootstrap
    seeds    = exp["bootstrap_seeds"]
    progress = DataFrame(name       = String[],
                         configfile = String[],
                         b = Int64[],
                         k = Int64[],
                         f = Array{Float64,1}[])
    for b in 1:length(seeds)
        info("Starting bootstrap $b/$(length(seeds)) for $configfile")
        rng = MersenneTwister(seeds[b])
        itrain, _   = Util.bootstrap_sample_indices(rng, length(y_train))
        idata,  oob = Util.bootstrap_sample_indices(rng, length(y_data))
        
        # true spectrum in this bootstrap sample (in-bag and oob)
        for (truth_name, indices) in zip([TRUESPEC_IB, TRUESPEC_OOB], [idata, oob])
            f_true = fit_pdf(y_data[indices], bins)
            push!(progress, [ truth_name, configfile, b, -1, f_true ])
        end
        
        # configure classifier and set up inspection
        tp = train_and_predict_proba(Util.classifier_from_config("conf/skl/nb.yml")) # from CherenkovDeconvolution.Sklearn
        inspect = (f, k, chi2s, alpha) -> push!(progress, [ name, configfile, b, k, f ])
                                                            
        # deconvolve the configured data set, representing one amount of smearing.
        dsea(X_data[idata, :], X_train[itrain, :], y_train[itrain], tp, bins,
             fixweighting = false, K = K, inspect = inspect)
    end
    return progress
    
end


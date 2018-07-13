"""
    weightfix(configfile="conf/job/weightfix.yml")

Experiments on the corrected re-weighting of training examples.
"""
function weightfix(configfile::String="conf/job/weightfix.yml")
    
    # read configuration
    c = load_file(configfile)
    B = c["num_bootstraps"]
    K = c["dsea_iterations"]
    outfile  = c["outfile"]
    datasets = c["datasets"]
    experiments = c["experiments"]
    
    # generate seeds for experiments and bootstrap
    srand(c["seed"])
    useeds = [ rand(UInt32) for _ in 1:maximum(map(d -> length(d["f_train"]), c["datasets"])) ]
    pushseeds!(experiments, B)
    
    # parallel execution
    info("Running experiments on $(length(datasets)) data sets, independently on $(nworkers()) worker(s)")
    df = vcat(pmap(dataset_conf -> _weightfix_job(dataset_conf, experiments, K, useeds), datasets)...)
    
    # output
    writetable(outfile, df)
    info("Results written to $outfile")
    df
    
end

# independent subroutine of weightfix
function _weightfix_job( dataset_conf :: Dict{Any,Any},
                         experiments  :: Array{Dict{Any,Any},1},
                         K            :: Int,
                         useeds       :: Array{UInt32,1} )
    # read data
    dataset_id = dataset_conf["id"]
    dataset = Data.dataset(dataset_id, nobs = -1, nobs_train = -1) # full data set
    discr   = Data.discretizer(dataset)
    bins    = Data.bins(discr)
    X_full  = Data.X_data(dataset)
    y_full  = encode(discr, Data.y_data(dataset))
    
    X_auxtrain, y_auxtrain = try
        (Data.X_train(dataset), encode(discr, Data.y_train(dataset)))
    catch (nothing, nothing) end
    
    # additional config
    skconfig      = joinpath("conf/skl", dataset_conf["skconfig"] * ".yml")
    f_train_confs = dataset_conf["f_train"] # string array
    
    # conduct experiments
    df = DataFrame()
    for (f_train_conf, useed) in zip(f_train_confs, useeds[1:length(f_train_confs)])
        X_data, y_data, X_train, y_train = Util.shuffle_split_subsample( X_full, y_full,
                                                                         X_auxtrain, y_auxtrain,
                                                                         f_train_conf,
                                                                         seed = useed )
        df = vcat(df, map(exp -> _weightfix_exp(exp, dataset_id, K, f_train_conf,
                                                X_data, y_data, X_train, y_train,
                                                bins, skconfig),  experiments)...)
    end
    return df
    
end

# single experiment
function _weightfix_exp{TN<:Number, TI<:Int}( exp          :: Dict{Any,Any}, 
                                              dataset_id   :: String,
                                              K            :: Int,
                                              f_train_conf :: String,
                                              X_data       :: Matrix{TN},
                                              y_data       :: Array{TI,1}, 
                                              X_train      :: Matrix{TN},
                                              y_train      :: Array{TI,1},
                                              bins         :: AbstractArray{TI,1},
                                              skconfig     :: String )
    
    # read basic config
    srand(exp["seed"])
    seeds        = exp["bootstrap_seeds"]
    fixweighting = exp["fixweighting"]
    name         = exp["name"]
    
    # bootstrap
    progress = DataFrame(name    = String[],
                         dataset = String[],
                         f_train = String[],
                         b       = Int64[],
                         k       = Int64[],
                         f       = Array{Float64, 1}[])
    for b in 1:length(seeds)
        info("Starting bootstrap $b/$(length(seeds)) of $name on $(dataset_id) data",
             " with $(f_train_conf) training set")
        rng = MersenneTwister(seeds[b])
        itrain, _   = Util.bootstrap_sample_indices(rng, length(y_train))
        idata,  oob = Util.bootstrap_sample_indices(rng, length(y_data))
        
        # true spectrum in this bootstrap sample (in-bag and oob)
        for (truth_name, indices) in zip([TRUESPEC_IB, TRUESPEC_OOB], [idata, oob])
            f_true = fit_pdf(y_data[indices], bins)
            push!(progress, [ truth_name, dataset_id, f_train_conf, b, -1, f_true ])
        end
        
        # configure classifier and set up inspection
        tp = train_and_predict_proba(Util.classifier_from_config(skconfig)) # from CherenkovDeconvolution.Sklearn
        inspect = (f, k, chi2s, alpha) -> push!(progress, [ name, dataset_id, f_train_conf, b, k, f ])
        
        # deconvolve the configured data set, representing one amount of smearing.
        dsea(X_data[idata, :], X_train[itrain, :], y_train[itrain], tp, bins,
             fixweighting = fixweighting, K = K, inspect = inspect)
    end
    return progress
    
end


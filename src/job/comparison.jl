"""
    comparison(configfile)

Comparative evaluation of RUN, IBU, and DSEA.
"""
function comparison(configfile::String)
    
    # read configuration
    c = load_file(configfile)
    B = c["num_bootstraps"]
    outfile     = c["outfile"]
    experiments = c["experiments"]
    
    # generate seeds for experiments and bootstraps
    srand(c["seed"])
    pushseeds!(experiments, B)
    for (i_exp, exp) in enumerate(experiments) # also add an info string
        exp["info"] = "$(i_exp)/$(length(experiments))"
    end
    
    # read data
    dataset  = Data.dataset(c["dataset"]["id"], nobs = -1, nobs_train = -1) # full data set
    discr    = Data.discretizer(dataset)
    X_full   = Data.X_data(dataset)
    y_full_c = Data.y_data(dataset) # continuous values, not indices!
    
    X_auxtrain, y_auxtrain_c = try
        (Data.X_train(dataset), Data.y_train(dataset))
    catch (nothing, nothing) end
    
    # more data-dependent configuration
    skconfig     = joinpath("conf/skl", c["dataset"]["skconfig"] * ".yml")
    f_train_conf = c["dataset"]["f_train"]
    
    # split into training and observed data sets without discretizing y_train
    X_data, y_data_c, X_train, y_train_c = Util.shuffle_split_subsample( X_full, y_full_c,
                                                                         X_auxtrain, y_auxtrain_c,
                                                                         f_train_conf,
                                                                         discretizer = discr )
    
    # parallel execution
    info("Starting $(length(experiments)) experiments on $(nworkers()) worker(s).")
    df = vcat(pmap(exp -> _comparison_job(exp, X_data, y_data_c, X_train, y_train_c, discr,
                                          f_train_conf, skconfig),
                   experiments)...)
    
    # output
    writetable(outfile, df)
    info("Results written to $outfile")
    df
    
end


# 
# Jobs of the different methods (DSEA, IBU, RUN) are delegated to the dedicated job
# functions implemented below. You can find the delegating function _comparison_job at the
# very end of this file.
# 


# progress DataFram of each experiment
_new_comparison_progress() = DataFrame( name           = String[],
                                        method         = String[], # deconvolution method
                                        fixweighting   = Bool[],   # weighting fix (DSEA)
                                        stepsize       = String[], # stepsize method (DSEA)
                                        smoothing      = String[], # smoothing method (IBU)
                                        discretization = String[], # discretization method (RUN)
                                        J              = Int64[],  # number of clusters (IBU, RUN)
                                        f_train        = String[],
                                        b              = Int64[],
                                        k              = Int64[],
                                        chi2s          = Float64[],
                                        f              = Array{Float64, 1}[] )


# independent sub-routine for DSEA
function _comparison_job_dsea{TN<:Number,TL<:Number}( exp       :: Dict{Any,Any},
                                                      X_data    :: Matrix{TN},
                                                      y_data_c  :: Array{TL,1}, 
                                                      X_train   :: Matrix{TN},
                                                      y_train_c :: Array{TL,1},
                                                      discr     :: AbstractDiscretizer,
                                                      f_train   :: String,
                                                      skconfig  :: String )
    # read basic configuration
    srand(exp["seed"])
    name            = exp["name"]
    seeds           = exp["bootstrap_seeds"]
    epsilon         = exp["epsilon"]
    K               = exp["K"]
    fixweighting    = exp["fixweighting"]
    stepsize_conf   = exp["stepsize"] # nested config
    stepsize_method = stepsize_conf["method"]
    
    # discretize the target quantity
    y_data  = encode(discr, y_data_c)
    y_train = encode(discr, y_train_c)
    bins    = Data.bins(discr)
    
    # bootstrap
    progress = _new_comparison_progress()
    for b in 1:length(seeds)
        info("Starting bootstrap $b/$(length(seeds)) of experiment $(exp["info"]) on $name")
        rng = MersenneTwister(seeds[b])
        itrain, _   = Util.bootstrap_sample_indices(rng, length(y_train))
        idata,  oob = Util.bootstrap_sample_indices(rng, length(y_data))
        
        # true spectrum in this bootstrap sample (in-bag and oob)
        for (truth_name, y) in zip([ TRUESPEC_IB, TRUESPEC_OOB, TRAINSPEC ], 
                                   [ y_data[idata], y_data[oob], y_train[itrain] ])
            push!(progress, [ truth_name, # name
                              "",         # method
                              false,      # fixweighting
                              "",         # stepsize
                              "",         # smoothing
                              "",         # discretization
                              -1,         # J
                              f_train,    # f_train
                              b,          # b
                              -1,         # k
                              NaN,        # chi2s
                              fit_pdf(y, bins) ]) # f
        end
        
        # line search strategy
        alpha = _stepsize_alpha(stepsize_conf, X_data[idata, :], X_train[itrain, :],
                                y_train[itrain], bins, zeros(0))
        
        # configure classifier and set up inspection
        tp = train_and_predict_proba(Util.classifier_from_config(skconfig)) # from CherenkovDeconvolution.Sklearn
        inspect = (f, k, chi2s, alpha) -> push!(progress, [ name, "dsea",
                                                            fixweighting,
                                                            stepsize_method,
                                                            "",  # smoothing
                                                            "",  # discretization
                                                            -1,  # J
                                                            f_train, b,
                                                            k, chi2s, f ])
        # start DSEA
        dsea(X_data[idata, :], X_train[itrain, :], y_train[itrain], tp, bins,
             alpha = alpha, fixweighting = fixweighting, epsilon = epsilon, K = K, inspect = inspect)
    end
    return progress
    
end


# independent sub-routine for IBU
function _comparison_job_ibu{TN<:Number,TL<:Number}( exp       :: Dict{Any,Any},
                                                     X_data    :: Matrix{TN},
                                                     y_data_c  :: Array{TL,1}, 
                                                     X_train   :: Matrix{TN},
                                                     y_train_c :: Array{TL,1},
                                                     discr     :: AbstractDiscretizer,
                                                     f_train   :: String,
                                                     skconfig  :: String )
    # read basic configuration
    srand(exp["seed"])
    name    = exp["name"]
    seeds   = exp["bootstrap_seeds"]
    epsilon = exp["epsilon"]
    K       = exp["K"]
    J       = exp["num_clusters"]
    smoothing_conf   = exp["smoothing"]
    smoothing_method = smoothing_conf["method"]
    
    # discretize the target quantity
    y_data  = encode(discr, y_data_c)
    y_train = encode(discr, y_train_c)
    bins    = Data.bins(discr)
    
    # bootstrap
    progress = _new_comparison_progress()
    for b in 1:length(seeds)
        info("Starting bootstrap $b/$(length(seeds)) of experiment $(exp["info"]) on $name")
        rng = MersenneTwister(seeds[b])
        itrain, _   = Util.bootstrap_sample_indices(rng, length(y_train))
        idata,  oob = Util.bootstrap_sample_indices(rng, length(y_data))
        
        # 
        # Note that the bootstrap truth is not written again, here. It is certain that
        # writing it would be redundant to what is already written by the DSEA job.
        # 
        
        # discretize the feature space with a decision tree from CherenkovDeconvolution.Sklearn
        discr_x = TreeDiscretizer(X_train[itrain, :], y_train[itrain], J)
        x_data  = encode(discr_x, X_data[idata,   :])
        x_train = encode(discr_x, X_train[itrain, :])
        
        # prepare smoothing function
        smoothing_func = Base.identity # assume smoothing_method == "none"
        if smoothing_method == "polynomial" # function from CherenkovDeconvolution.Util
            smoothing_func = polynomial_smoothing(smoothing_conf["order"], false)
        end
        
        # Iterative Bayesian Unfolding
        inspect = (f, k, chi2s) -> push!(progress, [ name, "ibu",
                                                     false, # fixweighting
                                                     "",    # stepsize method
                                                     smoothing_method,
                                                     "",    # discretization
                                                     J, f_train, b,
                                                     k, chi2s, f ])
        ibu(x_data, x_train, y_train[itrain], bins,
            smoothing = smoothing_func, epsilon = epsilon, K = K, inspect = inspect)
    end
    return progress
    
end


# independent sub-routine for RUN
function _comparison_job_run{TN<:Number,TL<:Number}( exp       :: Dict{Any,Any},
                                                     X_data    :: Matrix{TN},
                                                     y_data_c  :: Array{TL,1}, 
                                                     X_train   :: Matrix{TN},
                                                     y_train_c :: Array{TL,1},
                                                     discr     :: AbstractDiscretizer,
                                                     f_train   :: String,
                                                     skconfig  :: String )
    # read basic configuration
    srand(exp["seed"])
    name    = exp["name"]
    seeds   = exp["bootstrap_seeds"]
    epsilon = exp["epsilon"]
    K       = exp["K"]
    J       = exp["num_clusters"]
    f_dim   = length(Data.bins(discr)) # dimension of deconvolution results
    
    
    # discretize y_train, considering the expand/reduce configuration
    discr_conf   = exp["discretization"]
    discr_method = discr_conf["method"]
    discr_factor = discr_conf["factor"]
    if discr_method == "expand" && discr_factor > 1
        discr = expansion_discretizer(discr, discr_factor)
    elseif discr_method != "reduce" && discr_method != "expand"
        throw(ArgumentError("Discretization method '$(discr_method)' not supported"))
    end
    y_data  = encode(discr, y_data_c)
    y_train = encode(discr, y_train_c)
    bins    = Data.bins(discr)
    
    # choose n_df accordingly
    n_df = discr_method == "expand" ? f_dim : convert(Int, round(f_dim / discr_factor))
    
    
    # bootstrap
    progress = _new_comparison_progress()
    for b in 1:length(seeds)
        info("Starting bootstrap $b/$(length(seeds)) of experiment $(exp["info"]) on $name")
        rng = MersenneTwister(seeds[b])
        itrain, _   = Util.bootstrap_sample_indices(rng, length(y_train))
        idata,  oob = Util.bootstrap_sample_indices(rng, length(y_data))
        
        # 
        # Note that the bootstrap truth is not written again, here. It is certain that
        # writing it would be redundant to what is already written by the DSEA job.
        # 
        
        # discretize the feature space with a decision tree from CherenkovDeconvolution.Sklearn
        discr_x = TreeDiscretizer(X_train[itrain, :], y_train[itrain], J)
        x_data  = encode(discr_x, X_data[idata,   :])
        x_train = encode(discr_x, X_train[itrain, :])
        
        # prepare inspection with respect to expansion/reduction
        inspect = (f, k, ldiff, tau) -> push!(progress, [ name, "run",
                                                          false, # fixweighting
                                                          "",    # stepsize method
                                                          "",    # smoothing method
                                                          discr_method,
                                                          J, f_train, b,
                                                          k, ldiff, f ])
        if discr_factor > 1
            fun = inspect # copy pointer to inspection function
            if discr_method == "expand"
                inspect = inspect_expansion(fun, discr_factor)
            else
                inspect = inspect_reduction(fun, discr_factor)
            end
        end
        CherenkovDeconvolution.run(x_data, x_train, y_train[itrain], bins,
                                   n_df = n_df, epsilon = epsilon, K = K, inspect = inspect)
    end
    return progress
    
end


# delegate jobs to dedicated functions
_comparison_job(exp::Dict{Any,Any}, args...) = COMPARISON_JOB[exp["method"]](exp, args...)

COMPARISON_JOB = Dict( "dsea" => _comparison_job_dsea,
                       "ibu"  => _comparison_job_ibu,
                       "run"  => _comparison_job_run )



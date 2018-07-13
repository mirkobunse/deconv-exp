SMART_CONTROL_UNIFORM       = "uniform"
SMART_CONTROL_REDISTRICTING = "redistricting"
SMART_CONTROL_FTRAIN_U      = "ftrain uniform"
SMART_CONTROL_FTRAIN_R      = "ftrain redistricting"

"""
    smart_control(configfile="conf/job/smart_control.yml")

Experiments on smart control of simulations with DSEA.
"""
function smart_control(configfile::String="conf/job/smart_control.yml")
    
    # read configuration
    c = load_file(configfile)
    B = c["num_bootstraps"]
    dataset_id = c["dataset"]["id"]
    skconfig   = joinpath("conf/skl", c["dataset"]["skconfig"] * ".yml")
    outfile    = c["outfile"]
    smart_control = c["smart_control"]
    
    # generate bootstrap seeds
    srand(c["seed"])
    bootstrap_seeds = [ rand(UInt32) for _ in 1:B ]
    
    # generate training set sizes
    n_values = convert.(Int64, round.(logspace( log10(smart_control["min"]),
                                                log10(smart_control["max"]),
                                                smart_control["K"] )))
    
    # read data
    dataset = Data.dataset(dataset_id, nobs = -1)
    discr   = Data.discretizer(dataset)
    bins    = Data.bins(discr)
    X_full  = Data.X_data(dataset)
    y_full  = encode(discr, Data.y_data(dataset))
    
    X_auxtrain, y_auxtrain = try
        (Data.X_train(dataset), encode(discr, Data.y_train(dataset)))
    catch (nothing, nothing) end
    
    # split into training and observed data sets (do not subsample!)
    X_data, y_data, X_train, y_train = Util.shuffle_split_subsample( X_full, y_full,
                                                                     nothing, nothing,
                                                                     "appropriate",
                                                                     nobs_train = -1 )
    
    # add auxiliary data (if available) to the training set
    if X_auxtrain != nothing && y_auxtrain != nothing
        X_train = vcat(X_train, X_auxtrain)
        y_train = vcat(y_train, y_auxtrain)
    end
    
    # parallel execution
    info("Starting $(length(bootstrap_seeds)) bootstrap iterations on $(nworkers()) worker(s).")
    df = vcat(pmap(t -> _smart_control_job(t[1], t[2], X_data, y_data, X_train, y_train,
                                           bins, skconfig, n_values),
                   collect(enumerate(bootstrap_seeds)))...)
    
    # output
    writetable(outfile, df)
    info("Results written to $outfile")
    df
    
end

# independent sub-routine of smart_control
function _smart_control_job{TN<:Number,TI<:Int}( b          :: Int,
                                                 seed       :: UInt32,
                                                 X_data     :: Matrix{TN},
                                                 y_data     :: Array{TI,1},
                                                 X_train    :: Matrix{TN},
                                                 y_train    :: Array{TI,1},
                                                 bins       :: AbstractArray{TI,1},
                                                 skconfig   :: String,
                                                 n_values   :: Array{Int,1} )
    # bootstrap indices
    srand(seed)
    itrain, _   = Util.bootstrap_sample_indices(length(y_train))
    idata,  oob = Util.bootstrap_sample_indices(length(y_data))
    
    # true spectrum in this bootstrap sample (in-bag and oob)
    progress = DataFrame( name = String[],
                          b = Int64[],
                          n = Int64[],
                          f = Array{Float64, 1}[] )
    for (truth_name, indices) in zip([TRUESPEC_IB, TRUESPEC_OOB], [idata, oob])
        f_true = fit_pdf(y_data[indices], bins)
        push!(progress, [ truth_name, b, -1, f_true ])
    end
    
    # limit data to bootstrap
    X_data = X_data[idata, :]
    y_data = y_data[idata] # not used from here on
    X_train = X_train[itrain, :]
    y_train = y_train[itrain]
    
    # separate validation set from training pool
    X_val = X_train[1:10000, :]
    y_val = y_train[1:10000]
    X_train = X_train[10001:end, :]
    y_train = y_train[10001:end]
    
    # configure classifier
    c  = Util.classifier_from_config(skconfig)
    tp = train_and_predict_proba(c)
    
    # prepare storage for simulation control
    y_pred_1 = Int64[] # predictions in first iteration
    y_pred_2 = Int64[] # predictions in second iteration
    f_2      = Float64[] # estimate in second iteration
    
    
    # ======================================================================================
    # 
    # Uniform sampling as a base line
    # 
    utrain = Util.subsample_uniform_indices(y_train, shuffle = true)
    f_est  = Float64[] # will be interpreted as a uniform prior
    
    for n in n_values # increase the training set size
        info("Uniform sampling with $n training examples..")
        ntrain  = utrain[1:min(n, length(utrain))]
        f_train = fit_pdf(y_train[ntrain], bins)
        f_est   = dsea(X_data, X_train[ntrain, :], y_train[ntrain], tp, bins, f_0 = f_est)
        push!(progress, [ SMART_CONTROL_UNIFORM, b, n, f_est ])
        push!(progress, [ SMART_CONTROL_FTRAIN_U, b, n, f_train ])
        
        # prepare the redistricting sampling, which equals the uniform sampling in the first two rounds
        if n in n_values[1:2]
            push!(progress, [ SMART_CONTROL_REDISTRICTING, b, n, f_est ])
            push!(progress, [ SMART_CONTROL_FTRAIN_R, b, n, f_train ])
            if n == n_values[1]
                y_pred_1 = ScikitLearn.predict(c, X_val) # store prediction
            else
                y_pred_2 = ScikitLearn.predict(c, X_val)
                f_2      = f_est # also store estimate
            end
        end
    end
    
    
    # ======================================================================================
    # 
    # Simulation control with DSEA and redistriction weight (starting from iteration 2)
    # 
    ntrain = utrain[1:(n_values[2])]            # indices of the last training set
    ipool  = setdiff(1:length(y_train), ntrain) # pool of the remaining indices
    f_est  = f_2                                # the last estimate
    y_pred = y_pred_2                           # the last prediction
    w = _redistricting_weights(y_pred, y_pred_1, y_val, bins) # weights from predictions
    
    for n in n_values[3:end] # increase the training set size
        info("Smart control with $n training examples..")
        
        # sample the desired number of new examples in each bin from the pool
        f_train = fit_pdf(y_train[ntrain], bins)
        desired = (n - length(ntrain)) .* normalizepdf(normalizepdf(f_est + w) - f_train)
        ntrain = vcat(ntrain, _pop_from_pool!(ipool, y_train, desired))
        
        # deconvolve
        f_est = dsea(X_data, X_train[ntrain, :], y_train[ntrain], tp, bins, f_0 = f_est)
        push!(progress, [ SMART_CONTROL_REDISTRICTING, b, n, f_est ])
        
        # also store the training set density
        f_train = fit_pdf(y_train[ntrain], bins) # update density
        push!(progress, [ SMART_CONTROL_FTRAIN_R, b, n, f_train ])
        
        # weight by the number of redistricted examples in each bin
        y_prev = y_pred # previous prediction
        y_pred = ScikitLearn.predict(c, X_val)
        w = _redistricting_weights(y_pred, y_prev, y_val, bins)
    end
    return progress
end

# return indices of desired examples from y_train, which are removed from the pool
function _pop_from_pool!(ipool, y_train, desired)
    itrain = Int64[]
    for (i, n) in enumerate(convert.(Int64, round.(desired)))
        itrain_i = find(y_train[ipool] .== i) # find candidates in REMAINING PART of pool
        if length(itrain_i) >= n
            itrain_i = itrain_i[1:n]
        else
            warn("Desired $n new examples from bin $i, but only $(length(itrain_i)) are available!")
        end
        itrain = vcat(itrain, itrain_i)
        deleteat!(ipool, itrain_i)
    end
    return itrain
end

# weight by the number of redistricted examples in each bin
function _redistricting_weights( y_pred :: AbstractArray{Int64,1},
                                 y_prev :: AbstractArray{Int64,1},
                                 y_true :: AbstractArray{Int64,1},
                                 bins   :: AbstractArray{Int64,1} )
    n_r = map(bins) do b # count the examples for which the prediction changed
        i_b = y_true .== b # indices of examples with label b
        sum(y_pred[i_b] .!= y_prev[i_b])
    end
    return n_r ./ sum(n_r)
end


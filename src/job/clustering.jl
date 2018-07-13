"""
    clustering(configfile="conf/job/clustering.yml")

Experiments on the difficulty of classical deconvolution.
"""
function clustering(configfile::String="conf/job/clustering.yml")
    
    # read configuration
    c = load_file(configfile)
    B = c["num_bootstraps"]
    outfile  = c["outfile"]
    experiments = c["experiments"]
    
    # generate seeds for experiments and bootstrap
    srand(c["seed"])
    pushseeds!(experiments, B)
    
    # expand number of clusters
    expanded = Dict[]
    for exp in experiments
        num_clusters = c["num_clusters"]
        cluster_size = logspace(log10(num_clusters["min"]), log10(num_clusters["max"]), num_clusters["num_steps"])
        exp["num_clusters"] =  Int64.(round.(cluster_size)) # prepare for expansion
        push!(expanded, Yaml.expand(exp, "num_clusters")...)
    end
    
    # read data, discretizing the target variable
    dataset = Data.Gaussian()
    discr   = Data.discretizer(dataset)
    bins    = Data.bins(discr)
    X_data  = Data.X_data(dataset)
    X_train = Data.X_train(dataset)
    y_data  = encode(discr, Data.y_data(dataset))
    y_train = encode(discr, Data.y_train(dataset))
    
    # parallel execution
    info("Running $(length(expanded)) experiments independently on $(nworkers()) worker(s)")
    df = vcat(pmap(exp -> _clustering_job(exp, X_data, y_data, X_train, y_train, bins), expanded)...)
    
    # output
    writetable(outfile, df)
    info("Results written to $outfile")
    df
    
end

function _clustering_job( exp     :: Dict{Any,Any},
                          X_data  :: Matrix{Float64},
                          y_data  :: Array{Int64,1},
                          X_train :: Matrix{Float64},
                          y_train :: Array{Int64,1},
                          bins    :: AbstractArray{Int64,1} )
    srand(exp["seed"])
    
    # read basic configuration
    clustering = exp["clustering"]
    method     = clustering["method"]
    J          = exp["num_clusters"]
    seeds      = exp["bootstrap_seeds"]
    
    # meaningful name of strategy for plots
    name = if method == "univariate"
               "univariate"
           elseif method == "tree" && clustering["criterion"] == "gini"
               "Gini index"
           elseif method == "tree" && clustering["criterion"] == "entropy"
               "entropy"
           elseif method == "kmeans"
               "\$k\$-means"
           else  error("method '$method' unavailable")  end
    
    # bootstrap
    progress = DataFrame(name   = String[],
                         method = String[],
                         b      = Int64[],
                         J      = Int64[],
                         cond_R = Float64[],
                         f      = Array{Float64, 1}[])
    for b in 1:length(seeds)
        info("Starting bootstrap $b/$(length(seeds)) for strategy '$name' and J = $J")
        
        rng = MersenneTwister(seeds[b])
        itrain, _   = Util.bootstrap_sample_indices(rng, length(y_train))
        idata,  oob = Util.bootstrap_sample_indices(rng, length(y_data))
        
        # true spectrum in this bootstrap sample (in-bag and oob)
        for (truth_name, indices) in zip([TRUESPEC_IB, TRUESPEC_OOB], [idata, oob])
            f_true = fit_pdf(y_data[indices], bins)
            push!(progress, [ truth_name, "", b, -1, -1, f_true ])
        end
        
        # discretize observables
        bins_x  = Int64[]
        x_data  = Int64[]
        x_train = Int64[]
        if method == "univariate"
            # LinearDiscretizer from Discretizers.jl
            xmin, xmax = extrema(X_train[itrain, 1])
            discr_x = LinearDiscretizer(linspace(xmin, xmax, J+1))
            bins_x  = Data.bins(discr_x)
            x_data  = encode(discr_x, X_data[idata,   1])
            x_train = encode(discr_x, X_train[itrain, 1])
            
        elseif method == "tree"
            # TreeDiscretizer from CherenkovDeconvolution.Sklearn
            discr_x = TreeDiscretizer(X_train[itrain, :], y_train[itrain], J, clustering["criterion"])
            bins_x  = Sklearn.bins(discr_x)
            x_data  = encode(discr_x, X_data[idata,   :])
            x_train = encode(discr_x, X_train[itrain, :])
            
        elseif method == "kmeans"
            # KMeansDiscretizer from CherenkovDeconvolution.Sklearn
            discr_x = KMeansDiscretizer(X_train[itrain, :], J)
            bins_x  = Sklearn.bins(discr_x)
            x_data  = encode(discr_x, X_data[idata,   :])
            x_train = encode(discr_x, X_train[itrain, :])
            
        else
            error("Clustering method '$method' not available")
        end
        
        # estimate transfer and observed distribution (in-bag)
        R = fit_R(y_train[itrain], x_train, bins_x=bins_x, bins_y=bins)
        g = fit_pdf(x_data, bins_x)
        
        # minimum-norm least-squares fit = Moore-Penrose pseudoinverse solution
        f = pinv(R) * g
        
        cond_R = try  cond(R)  catch err  warn(err); Inf  end
        push!(progress, [ name, method, b, J, cond_R, f ])
        
    end
    return progress
    
end

"""
    gridsearch(configfile="conf/job/gridsearch.yml")

Grid search experiment finding suitable parameters for a random forest on FACT data.
"""
function gridsearch(configfile::String="conf/job/gridsearch.yml")
    
    # configuration
    c = load_file(configfile)
    srand(c["seed"])
    dataset  = Data.dataset(c["dataset"])
    skconfig = joinpath("conf/skl", c["skconfig"] * ".yml")
    outfile  = c["outfile"]
    parameters   = Dict(zip(map(Symbol, keys(c["parameters"])),
                          map(x -> eval(parse(x)), values(c["parameters"]))))
    
    # subsample uniform data set
    discr = Data.discretizer(dataset)
    X, y  = Util.subsample_uniform(Data.X_data(dataset),
                                   encode(discr, Data.y_data(dataset)))
    
    # perform grid search
    bestparams = gridsearch(X, y, skconfig, parameters)
    
    # update parameters in skconfig
    skconfig = load_file(skconfig)
    for k in keys(bestparams)
        skconfig["parameters"][string(k)] = string(bestparams[k]) * " # " * c["parameters"][string(k)]
    end
    
    # write skconfig to file
    info("Writing configuration to $outfile")
    prefix = """
	# 
	# $outfile
	# 
	# Automatically generated with Job.gridsearch(\"$configfile\").
	# 
	"""
    write_file(outfile, skconfig, prefix)
    info("Results written to $outfile")
    
end


"""
    gridsearch(X, y, configfile, params)

Perform a grid search on the feature matrix `X` and labels `y`, applying the classifier
configured in the `configfile`. `params` is a dictionary mapping from parameter symbols to
value arrays.
"""
function gridsearch{TN<:Number, TI<:Int}(X::Matrix{TN}, y::Array{TI,1}, configfile::String,
                                         params::Dict)
    
    num_combinations = (*(map(length, values(params))...))
    _reset_gridsearch_progress(3 * num_combinations) # 3 folds is the default
    info("Performing a 3-fold CV grid search on $num_combinations parameter combinations.")
    
    # perform grid search
    classifier = Util.classifier_from_config(configfile)
    search = ScikitLearn.GridSearch.GridSearchCV(classifier, params,
                                                 scoring = _gridsearch_scoring)
    ScikitLearn.fit!(search, X, y)
    
    info("Best score is $(search.best_score_)")
    return search.best_params_
    
end

_gridsearch_progress = [ 0, 0 ] # total and current number

# reset progress to i total iterations
function _reset_gridsearch_progress(i::Int)
    deleteat!(_gridsearch_progress, [1, 2])
    push!(_gridsearch_progress, i)
    push!(_gridsearch_progress, 1)
end

# advance and print last iteration
function _advance_gridsearch_progress()
    i = pop!(_gridsearch_progress)
    push!(_gridsearch_progress, i + 1)
    return "$i/$(_gridsearch_progress[1])"
end

# score a parameter set - high score is good, so it is the negative error
function _gridsearch_scoring{TN<:Number, TI<:Int}(clf::PyObject,
                                                  X::Matrix{TN},
                                                  y_true::Array{TI,1};
                                                  sample_weight=nothing)
    info("Grid search evaluation ", _advance_gridsearch_progress())
    y_pred = ScikitLearn.predict(clf, X)
    return - weighted_sum(abs.(y_true - y_pred), sample_weight; normalize=true)
end


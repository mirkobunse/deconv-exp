module Job


using DataFrames, Optim, MLDataUtils, Distances, Distributions, Discretizers
using ScikitLearn, ScikitLearnBase.weighted_sum, PyCall

using ComfyCommons: info, warn, Git, Yaml
using CherenkovDeconvolution
using CherenkovDeconvolution: Util, Sklearn

import Util, Data


# names for true spectra
TRUESPEC_IB  = "TRUE_SPECTRUM_IN_BAG"
TRUESPEC_OOB = "TRUE_SPECTRUM_OOB"
TRAINSPEC    = "TRAIN_SPECTRUM"


"""
    run(configfile)

Hand the `configfile` to the Job method of which the name is configured by the property
`job` in the `configfile`.
"""
function run(configfile::String)
    c = load_file(configfile)
    funname = "Job." * c["job"]
    
    info("Calling $funname(\"$configfile\")")
    fun = eval(parse(funname))
    fun(configfile) # function call
end


"""
    pushseeds!(experiments, B)

Add a random experiment seed and `B` seeds for the bootstrap samples to all configurations
in `experiments`. The bootstrap samples are equal in all experiments for comparability.

The experiment seed is stored in the property `seed`, the bootstrap seeds are an array in
the `bootstrap_seeds` property. When the global random number generator is seeded with
`srand()`, the result of this `pushseeds!()` is deterministic.
"""
function pushseeds!(experiments::AbstractArray{Dict{Any, Any}, 1}, B::Int)
    bootstrapseeds = [ rand(UInt32) for _ in 1:B ] # equal in all experiments
    for exp in experiments
        exp["seed"] = rand(UInt32) # individual seed
        exp["bootstrap_seeds"] = bootstrapseeds
    end
end


# job implementations
include("job/smearing.jl")   # different amounts of smearing in data
include("job/clustering.jl") # the difficulty of classical deconvolution
include("job/gridsearch.jl") # suitable parameters of a random forest on FACT data
include("job/weightfix.jl")  # corrected re-weighting of training examples
include("job/stepsize.jl")   # step size extension
include("job/smoothing.jl")  # smoothing extension
include("job/expand_reduce.jl") # expansion/reduction in DSEA
include("job/comparison.jl")    # comparative evaluation of RUN, IBU, and DSEA
include("job/time_series_contributions.jl") # contributions for time series analyses
include("job/smart_control.jl") # smart control of simulations with DSEA


end


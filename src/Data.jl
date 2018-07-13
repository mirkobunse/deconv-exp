module Data


using DataFrames, Distributions, Discretizers, HDF5
using ComfyCommons: info, warn, Yaml

import Util


DISCRETIZATION_Y = "discretization_y" # property in config files


"""
    LinearDiscretizer(configfile, key=""; kwargs...)

Configure a `LinearDiscretizer` object with the YAML `configfile`.

The YAML file should specify the properties `min`, `max`, and `num_bins`. If `key` is not
empty, these properties are looked for under the given `key`. All of these properties can be
substituted with the keyword arguments, e.g. with `min = 0.4`.

Example configuration file
--------------------------

    # myconf.yml
    mykey:
      min: 0.3
      max: 0.8
      num_bins: 3
"""
function Discretizers.LinearDiscretizer(configfile::String, key::String = DISCRETIZATION_Y;
                                        kwargs...)
    # load config
    c = load_file(configfile)
    if length(key) > 0
        c = c[key]
    end
    
    # merge arguments into config
    for (k, v) in kwargs
        c[string(k)] = v
    end
    
    return LinearDiscretizer(linspace(c["min"], c["max"], c["num_bins"]+1))
end

"""
    bins(ld)

Obtain the bin indices of the LinearDiscretizer (or DataSet object) `ld`.
"""
bins(ld::LinearDiscretizer) = indices(bincenters(ld), 1) # DataSet impl is below


# DataSet type implementations
abstract type DataSet end   # abstract supertype
include("data/gaussian.jl") # toy data
include("data/fact.jl")     # FACT telescope data
include("data/magic.jl")    # MAGIC telescope data

"""
    X_data(d)

Return the feature matrix `X` of observed data in the DataSet `d`.
"""
X_data(d::DataSet) = throw(ArgumentError("X_data is not implemented for $(typeof(d))"))

"""
    y_data(d)

Return the target value array `y` of observed data in the DataSet `d`.
"""
y_data(d::DataSet) = throw(ArgumentError("y_data is not implemented for $(typeof(d))"))

"""
    X_train(d)

Return the feature matrix `X` of training data in the DataSet `d`.
"""
X_train(d::DataSet) = throw(ArgumentError("X_train is not implemented for $(typeof(d))"))

"""
    y_train(d)

Return the target value array `y` of training data in the DataSet `d`.
"""
y_train(d::DataSet) = throw(ArgumentError("y_train is not implemented for $(typeof(d))"))

"""
    discretizer(d)

Return the target value discretizer for the DataSet `d`.
"""
discretizer(d::DataSet) = throw(ArgumentError("discretizer is not implemented for $(typeof(d))"))

bins(d::DataSet) = bins(discretizer(d))

"""
    dataset(id)

Return the DataSet object with the given `id`.
"""
dataset(id::String, args...; kwargs...) = DATASET_TYPES[id](args...; kwargs...)

DATASET_TYPES = Dict( "gaussian" => Gaussian,
                      "fact"     => Fact,
                      "magic"    => Magic )


end


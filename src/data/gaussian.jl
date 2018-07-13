GAUSSIAN_GROUP_DATA  = "data"
GAUSSIAN_GROUP_TRAIN = "train"

"""
    type Gaussian <: DataSet
    
    Gaussian([configfile = "conf/data/gaussian.yml";
              nobs_data = 50000,
              nobs_train = 100000,
              readdata = true])

Toy data set of which the target distribution is a gaussian mixture model with two
components.
"""
type Gaussian <: DataSet
    
    configfile::String
    X_data::Matrix{Float64}
    y_data::Array{Float64, 1}
    X_train::Matrix{Float64}
    y_train::Array{Float64, 1}
    discretizer::LinearDiscretizer
    
    function Gaussian(configfile::String="conf/data/gaussian.yml";
                      nobs::Int = 50000,
                      nobs_train::Int = 100000,
                      readdata::Bool = true)
        # read data
        X_data,  y_data  = _Xy_gaussian(configfile, GAUSSIAN_GROUP_DATA,  nobs, readdata)
        X_train, y_train = _Xy_gaussian(configfile, GAUSSIAN_GROUP_TRAIN, nobs_train, readdata)
        return new(configfile, X_data, y_data, X_train, y_train,
                   LinearDiscretizer(configfile))
    end
    
end

# implementation of interface
X_data(d::Gaussian) = d.X_data
y_data(d::Gaussian) = d.y_data
X_train(d::Gaussian) = d.X_train
y_train(d::Gaussian) = d.y_train
discretizer(d::Gaussian) = d.discretizer

# read gaussian data
_Xy_gaussian(configfile::String, group::String, nobs::Int, readdata::Bool) =
    if readdata
        h5open(load_file(configfile)["datafile"], "r") do file
            group_y = file[group * "/y"]
            y = nobs > 0 ? group_y[1:nobs] : read(group_y)
            X = hcat(map(j -> begin
                                  group_x = file[group * "/x" * string(j)]
                                  nobs > 0 ? group_x[1:nobs] : read(group_x)
                              end,
                         1:(length(file[group])-1))...)
            X, y # return value
        end
    else
        zeros(0, 0), zeros(0) # empty X and y
    end


"""
    generate_gaussian(configfile = "conf/data/gaussian.yml")

Generate data as specified in the `configfile`. The generation strategy strictly follows the
implementation provided with [ruhe2013data], which is different from the description given
in that reference. The parameter choice is more convent here.
"""
function generate_gaussian(configfile = "conf/data/gaussian.yml")
    c = load_file(configfile)
    
    # set up parameters
    ca = c["attributes"]
    srand(ca["seed"])
    p = DataFrame(attribute = map(j -> Symbol("x" * string(j)), 1:ca["num_attributes"]),
                  num = collect(1:ca["num_attributes"]))
    for param in ["offset", "factor", "gamma", "sigma"]
        p[Symbol(param)] = rand(ca["num_attributes"]) .* (ca[param]["max"] - ca[param]["min"]) .+ ca[param]["min"]
    end
    
    # generate true quantity
    cy = c["y"]
    min_y = c[DISCRETIZATION_Y]["min"]
    max_y = c[DISCRETIZATION_Y]["max"]
    ymodel = Truncated(MixtureModel(Normal, [(cy["mean1"], cy["sigma1"]),
                                             (cy["mean2"], cy["sigma2"])],
                                    [cy["weight1"], cy["weight2"]]), cy["min"], cy["max"])
    srand(c["event_seed"])
    data  = DataFrame(y = rand(ymodel, c["num_events"])) # gaussian mixture
    train = DataFrame(y = rand(c["num_events"]) .* (max_y - min_y) .+ min_y) # uniform
    println("Generated true quantity y...")
    
    # 
    # generate secondary attributes
    # 
    # description in [ruhe2013data]: x_{n, i} = N(offset_i + factor_i * y_n ^ gamma_i, sigma_i)
    # implementation (here & there): x_{n, i} = offset_i + factor_i * N(y_n, sigma_i) ^ gamma_i
    #
    for j in 1:ca["num_attributes"]
        for df in [data, train]
            smeared = (p[j,:sigma] > 0) ? [ rand(Normal(y, p[j,:sigma])) for y in df[:y] ] : df[:y]
            df[Symbol("x_"*string(j))] = ((smeared .^ p[j,:gamma]) .* p[j,:factor]) .+ p[j,:offset]
        end
        println("Generated attribute $j/$(ca["num_attributes"])...")
    end
    
    # sort secondary attributes based on their correlation with true quantity
    cors = map(j -> cor(data[:y], data[Symbol("x_" * string(j))]), 1:ca["num_attributes"])
    attperm = sortperm(cors, rev=true) # permutation that orders attributes by their correlation with x
    invperm = zeros(attperm)
    for i=1:length(invperm)
        invperm[attperm[i]] = i
    end
    rename!(data,  setdiff(names(data),  [:y]), map(j -> Symbol("x" * string(j)), invperm))
    rename!(train, setdiff(names(train), [:y]), map(j -> Symbol("x" * string(j)), invperm))
    p[:attribute] = p[invperm, :attribute]
    p[:num] = p[invperm, :num]
    sort!(p, cols=[:num])
    p = p[:, setdiff(names(p), [:num])]
    
    # write to HDF5 file
    cols = vcat([:y], p[:attribute])
    h5open(c["datafile"], "w") do file
        for col in cols
            write(file, "data/"*string(col), convert(Array{Float64,1}, data[col]))
            write(file, "train/"*string(col), convert(Array{Float64,1}, train[col]))
        end
        for param in names(p)
            write(file, "parameters/"*string(param),
                  param == :attribute ? convert(Array{String,1}, map(att -> string(att), p[param]))
                                      : convert(Array{Float64,1}, p[param]))
        end
    end
    
end

"""
    parameters_gaussian(configfile = "conf/data/gaussian.yml")

Read parameters of previously generated data.
"""
function parameters_gaussian(configfile::String = "conf/data/gaussian.yml")
    df = DataFrame()
    h5open(load_file(configfile)["datafile"], "r") do file
        for param in ["attribute", "offset", "factor", "gamma", "sigma"]
            column = read(file["parameters/"*param])
            df[Symbol(param)] = param == "attribute" ? convert(Array{Symbol,1}, map(att -> Symbol(att), column)) : column
        end
    end
    df
end


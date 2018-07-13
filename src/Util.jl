module Util


using DataFrames, MLDataUtils, Discretizers, Distances, ScikitLearn
using ComfyCommons: info, warn, Yaml

@sk_import naive_bayes : GaussianNB
@sk_import ensemble    : RandomForestClassifier
@sk_import calibration : CalibratedClassifierCV # for future work


DISTANCE_EPSILON = 1e-9 # min value of pdfs assumed for distance computations


"""
    classifier_from_config(configfile)

Obtain the result of the `classifier` function from a YAML configuration instead of
using function arguments.
"""
function classifier_from_config(configfile::AbstractString)
    c = load_file(configfile) # read config
    classname = c["classifier"]
    params    = get(c, "parameters", nothing) != nothing ? c["parameters"] : Dict{Symbol, Any}()
    calibrate = get(c, "calibrate",  false)
    return classifier(classname, calibrate;
                      zip(map(Symbol, keys(params)), values(params))...)
end

"""
    classifier(classname, calibrate = false; kwargs...)

Obtain a classifier to be used in `train_and_predict_proba`. The following values of
`classname` are available:

- GaussianNB
- RandomForestClassifier

The keyword arguments configure the corresponding class (see official scikit-learn doc).
"""
function classifier(classname::AbstractString,
                    calibrate::Bool = false;
                    kwargs...)
    Classifier = eval(parse(classname)) # constructor method
    classifier = Classifier(; kwargs...)
    if calibrate
        classifier = CalibratedClassifierCV(classifier, method="isotonic")
    end
    return classifier
end


"""
    shuffle_split_subsample(X_full, y_full, X_auxtrain, y_auxtrain, f_train[; seed,
                            nobs = 50000, nobs_train = 100000)

Shuffle the `X_full` and `y_full`, split this data into training and observed data sets, and
subsample the training part. Three values for `f_train` configure the kind of subsampling:

- `"appropriate"` means no subsampling and it thus maintains the inherent distribution of
  `X_full` and `y_full`.
- `"uniform"` undersamples a training set that is uniformly distributed in the target
  variable.
- `"auxiliary"` replaces the training part obtained by `X_full` and `y_full` with a shuffled
  subset of the auxiliary data provided with `X_auxtrain` and `y_auxtrain`.

**Returns**

- `X_data` the observed feature matrix
- `y_data` the target values of the observed data
- `X_train` the training feature matrix
- `y_train` the training labels

`nobs` and `nobs_train` optionally specify the maximum number of observations in these data
sets.
"""
function shuffle_split_subsample{TN<:Number, TL<:Number}(
            X_full      :: Matrix{TN},
            y_full      :: Array{TL,1},
            X_auxtrain  :: Union{Void,Matrix{TN}},
            y_auxtrain  :: Union{Void,Array{TL,1}},
            f_train     :: String;
            seed        :: Integer = convert(Int, rand(UInt32)),
            nobs        :: Int = 50000,
            nobs_train  :: Int = 100000,
            discretizer :: Union{Void,AbstractDiscretizer} = nothing )
    
    # shuffle and split
    urng   = MersenneTwister(seed)
    i_rand = randperm(urng, length(y_full))  # shuffled indices
    X_data = X_full[i_rand[1:nobs], :]        # first 50000 examples, by default
    y_data = y_full[i_rand[1:nobs]]
    X_train = X_full[i_rand[(nobs+1):end], :] # assume f_train=="appropriate" (original distribution)
    y_train = y_full[i_rand[(nobs+1):end]]
    
    # subsample or replace by auxiliary training data
    if f_train == "uniform"
        if discretizer == nothing
            X_train, y_train = subsample_uniform(X_train, y_train) # discrete y_train
        else
            X_train, y_train = subsample_uniform(X_train, y_train, discretizer) # continuous y_train
        end
    elseif f_train == "auxiliary"
        i_aux   = randperm(urng, length(y_auxtrain)) # shuffle auxiliary data, too
        X_train = X_auxtrain[i_aux, :]
        y_train = y_auxtrain[i_aux]
    elseif f_train != "appropriate"
        throw(ArgumentException("'$(f_train)' is not a legal value for f_train"))
    end
    
    # by default, limit the training data to 100000 examples
    if nobs_train > 0
        nobs_train = min(length(y_train), nobs_train)
        X_train = X_train[1:nobs_train, :]
        y_train = y_train[1:nobs_train]
    end
    
    return X_data, y_data, X_train, y_train
end

"""
    subsample_uniform(X, y[, d; shuffle = false])

Undersample a uniformly distributed data set consisting of the features matrix `X` and the
label array `y`.

If an AbstractDiscretizer `d` is specified, discretize the continuous values from `y` to
find its label while maintaining the continuous values of `y` in the result.
"""
function subsample_uniform{TN<:Number, TI<:Int}(X::Matrix{TN}, y::Array{TI,1}; shuffle::Bool=false)
    X, y = undersample((X', y), shuffle = shuffle)
    return X', convert(Array, y) # repair output of MLDataUtils.undersample
end

function subsample_uniform{TN<:Number, TL<:Number}(X :: Matrix{TN},
                                                   y :: Array{TL,1},
                                                   d :: AbstractDiscretizer;
                                                   shuffle::Bool=false)
    X, y = undersample(v -> encode(d, v), (X', y), shuffle = shuffle)
    return X', convert(Array, y) # repair output of MLDataUtils.undersample
end

"""
    subsample_uniform_indices(y[; shuffle = false])

Undersample a uniformly distributed data set with the label array `y`. Return the indices of
the sub-sample instead of the sub-sample itself.
"""
function subsample_uniform_indices{TI<:Int}(y::Array{TI,1}; shuffle::Bool=false)
    y = undersample(y, shuffle = shuffle)
    return y.indexes[1]
end


"""
    mdpa(a, b)

Minimum Distance of Pair Assignments (MDPA) [cha2002measuring] for ordinal pdfs `a` and `b`.
The MDPA is a special case of the Earth Mover's Distance [rubner1998metric] that can be
computed efficiently.
"""
function mdpa{T<:Number}(a::AbstractArray{T,1}, b::AbstractArray{T,1})
    if length(a) != length(b)
        throw("histograms have to have the same length")
    elseif !isapprox(sum(a), sum(b))
        throw("histograms have to have the same mass (difference is $(sum(a)-sum(b))")
    end
    
    # algorithm 1 in [cha2002measuring]
    prefixsum = 0.0
    distance  = 0.0
    for i in 1:length(a)
        prefixsum += a[i] - b[i]
        distance  += abs(prefixsum)
    end
    return distance / sum(a) # the normalization is a fix to the original MDPA
end

pairwise_mdpa{T<:Number}(a::AbstractArray{T,2}) = pairwise_mdpa(a, a)
function pairwise_mdpa{T<:Number}(a::AbstractArray{T,2}, b::AbstractArray{T,2})
    mat = zeros(size(a, 2), size(b, 2))
    for i in 1:size(a, 2), j in i:size(b, 2)
        mat[i,j] = mdpa(a[:,i], b[:,j])
    end
    return full(Symmetric(mat))
end

"""
    chi2p(a, b)

Pearson's Chi Square distance between the pdfs `a` and `b`, where `b` is the truth.
"""
chi2p{T<:Number}(a::AbstractArray{T,1}, b::AbstractArray{T,1}) =
        Distances.wsqeuclidean(a, b, 1 ./ max.(b, DISTANCE_EPSILON))

"""
    kl(a, b)

Kullback-Leibler divergence between the pdfs `a` and `b`, where `b` is the truth.
"""
kl{T<:Number}(a::AbstractArray{T,1}, b::AbstractArray{T,1}) =
        Distances.kl_divergence(a, max.(b, DISTANCE_EPSILON)) # KL divergence only defined for b_i != 0

"""
    chi2s(a, b)

Symmetric Chi Square distance between the pdfs `a` and `b`.
"""
chi2s{T<:Number}(a::AbstractArray{T,1}, b::AbstractArray{T,1}) = 2 * Distances.chisq_dist(a, b)


"""
    bootstrap_sample_indices([rng, ] N)

Return the training set and test set indices for a data pool of size `N`.
"""
bootstrap_sample_indices(N::Int64) = bootstrap_sample_indices(Base.GLOBAL_RNG, N)
function bootstrap_sample_indices(rng::AbstractRNG, N::Int64)
    i_train = rand(rng, 1:N, N)
    return i_train, setdiff(1:N, i_train) # tuple of training and test set indices
end


"""
    latex_e(n[, digits = 2, dollars = true])

Format the number `n` in scientific LaTeX format.
"""
latex_e(n::Float64, digits::Int=2; dollars::Bool=true) =
    if n != 0
        m = match(r"(\d\.\d+)e(-?\d+)", @sprintf("%e", n))
        b = round(parse(Float64, m[1]), digits)
        s = (b == 1.0 ? "" : "$b \\cdot ") * "10^{$(m[2])}"
        return dollars ? "\$$s\$" : s
    else
        return "0"
    end


"""
    binarytransfer(df, y; normalize=true)

Empirically estimate the transfer matrix from the `y` column in the DataFrame `df` to the
item index. The resulting normalized matrix `R` can be used to deconvolve in a next-neighbor
fashion.

**Currently unused function for future work**
"""
function binarytransfer{T <: Number}(df::DataFrame, y::Symbol;
                                     ylevels::AbstractArray{T, 1} = sort(unique(df[y])),
                                     normalize::Bool = true) # TODO change interface to Array types
    # construct binary matrix
    onehot = y_val::Float64 -> _diracimpulse(findfirst(ylevels, y_val), length(ylevels))
    R = reshape(vcat([ onehot(y_val) for y_val in df[y] ]...), (length(ylevels), size(df,1)))'
    
    if normalize
        return normalizetransfer!(R)
    else
        return R
    end
    
end

"""
    nextneighbor_pdf(data, train, features = names(df); batchsize = 100)

A histogram of next neighbors to be used together with the matrix returned by
`Util.binarytransfer()`.

**Caution:** Throws an `OutOfMemoryError` when the batch size is too large.

**Currently unused function for future work**
"""
function nextneighbor_pdf(data::DataFrame, train::DataFrame,
                          features::Array{Symbol,1} = names(df);
                          batchsize::Int64 = 100) # TODO change interface to Array types
    data  = data[:,  features]
    train = train[:, features]
    g = zeros(Int64, size(train, 1))
    for batch in batchview(data, size = batchsize)
        for i in _findnextneighbors(batch, train)
            g[i] += 1
        end
    end
    return g
end

_findnextneighbors(data::AbstractDataFrame, train::AbstractDataFrame) =
    mapslices(col -> findmin(col)[2],
              pairwise(Euclidean(), convert(Array, data)', convert(Array, train)'), 2)

# Return a vector of length `l` that is `1.0` at position `p` and `0.0` in all other dimensions.
_diracimpulse(p::Int, l::Int) = (d = zeros(l);  d[p] = 1.0;  d)

end

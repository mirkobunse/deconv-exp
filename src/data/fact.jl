"""
    type Fact <: DataSet
    
    Fact([configfile = "conf/data/fact.yml"; readdata = true, nobs = -1])

FACT telescope data set.

**Caution:** This data set is not available in public.
"""
type Fact <: DataSet
    
    configfile::String
    X_data::Matrix{Float64}
    y_data::Array{Float64, 1}
    discretizer::LinearDiscretizer
    
    function Fact(configfile::String="conf/data/fact.yml"; nobs::Int = -1,
                  readdata::Bool = true, kwargs...)
        
        # 
        # CAUTION: This data set is not available in public
        # 
        if readdata # else, everything is fine because nothing is read
            error("The data configured by $configfile is not available in public")
        end
        
        # 
        # If it was, the following code would be executed...
        # 
        
        if length(kwargs) > 0 # just warn, do not throw an error
            warn("DataSet configured by $configfile has no keyword argument $(kwargs[1][1])")
        end
        X_data, y_data = _Xy_fact(configfile, nobs, readdata)
        return new(configfile, X_data, y_data, LinearDiscretizer(configfile))
    end
    
end

# implementation of interface
X_data(d::Fact) = d.X_data
y_data(d::Fact) = d.y_data
discretizer(d::Fact) = d.discretizer

# read FACT data
_Xy_fact(configfile::String, nobs::Int, readdata::Bool) =
    if readdata
        df = readtable(load_file(configfile)["datafile"], nrows = nobs)
        y = Float64.(convert(Array, df[FACT_TARGET]))
        X = Float64.(convert(Array, df[:, setdiff(names(df), [FACT_TARGET])]))
        X, y # return value
    else
        zeros(0, 0), zeros(0) # empty X and y
    end

FACT_TARGET = :log10_energy


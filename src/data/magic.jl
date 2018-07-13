"""
    Magic([configfile = "conf/data/magic.yml"; nobs = -1])

Magic telescope data set.

**Caution:** This data set is not available in public.
"""
Magic(; kwargs...) = Fact("conf/data/magic.yml"; kwargs...)


"""
    fix_magic(configfile = "conf/data/magic.yml")

Fix the files obtained from docker-truee, where the original MAGIC .root files are converted
to multiple CSV files. This method renames the columns conveniently and adds the zenith
angle group to the attributes. Moreover, the logarithm of some attributes is computed.
Finally, a single CSV file is written to disk.
"""
function fix_magic(configfile::String = "conf/data/magic.yml")
    
    # read config
    c = load_file(configfile)
    renamedict = Dict(map(tup -> (Symbol(tup[1]), Symbol(tup[2])), [ c["rename"]... ]))
    df = DataFrame()
    
    # read all input files
    regex = Regex(basename(c["infiles"]))
    files = filter(p -> ismatch(regex, p), readdir(dirname(c["infiles"])))
    for file in files
        filepath = joinpath(dirname(c["infiles"]), file)
        info("Reading $filepath")
        za = match(regex, file)[1] # zenith angle group
        sdf = readtable(filepath, normalizenames=false)
        rename!(sdf, renamedict)
        df = vcat(df, sdf)
    end
    
    # apply logarithm
    for col in map(Symbol, c["logarithm"])
        newcol = Symbol("log10_$col")
        info("Generating $newcol")
        df[col] = log10(df[col])
        rename!(df, col, newcol)
    end
    
    info("Writing to $(c["datafile"])")
    writetable(c["datafile"], df) # write to CSV (HDF5 would not save much space)
end


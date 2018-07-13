"""
    clustering(metricsfile="res/metrics/clustering.csv")


"""
function clustering(metricsfile::String="res/metrics/clustering.csv"; full::Bool=false)
    initialize("clustering", metricsfile) # initialize PGFPlots
    df = readtable(metricsfile) # read metrics
    
    # select relevant strategies
    if !full
        df = df[df[:name] .!= "entropy", :]
        df[df[:name] .== "Gini index", :name] = "{\\sc Cart} tree" # rename strategy
    end
    
    # prepare ticks in plots
    xticks     = sort(unique(df[:J]))[unique(vcat(1:5:end, end))]
    xticklabel = "{\\pgfkeys{/pgf/fpu=true}\\pgfmathparse{exp(\\tick)}\\pgfmathprintnumber[fixed relative, precision=3]{\\pgfmathresult}\\pgfkeys{/pgf/fpu=false}}"
    
    # plot
    for ref in [ "oob", "ib" ], mkey in vcat(:cond_R, keys(METRICS)...)
        
        # metric column from reference string and metric (:cond_R without reference)
        metric = Symbol("$(ref)_$(mkey)") # default
        iscond = mkey == :cond_R         # bool is reused later
        if iscond
            if ref == "ib" # cond_R equal for ib and oob
                continue
            end
            metric = mkey   # key without reference
        end
        
        # aggregate mean and quantiles (5% and 95%)
        agg = aggregate_bootstrap(df, [:name, :J], metric)
        
        # file paths
        outfile_pdf = pdfpath(metricsfile, "_" * string(metric))
        outfile_tex = texpath(metricsfile, "_" * string(metric))
        info("Plotting to $(outfile_pdf) and $(outfile_tex)")
        
        # plots differ in y axis label
        ylabel = iscond ? "\$\\kappa_\\mathbf{R}\$" : METRIC_NAMES[mkey] * " to \$\\mathbf{f}\$"
        
        # axis with one layer per clustering
        plot = Axis( xlabel = "Number of clusters \$J\$",
                     ylabel = ylabel,
                     style  = join([ "progress axis",
                                     "scale = .9",
                                     "ymode = log",
                                     "xmode = log",
                                     "xtick = {$(join(xticks, ", "))}",
                                     "xticklabel = $xticklabel",
                                     "legend pos = north east" ], ", ") )
        for (i, name) in enumerate(unique(df[:name]))
            sdf = agg[agg[:name] .== name, :]
            push!(plot, Plots.Linear( sdf[:J], sdf[:y],
                                      errorBars   = ErrorBars( yplus  = sdf[:y_plus],
                                                               yminus = sdf[:y_minus]),
                                      style       = progress_style(i),
                                      legendentry = name ))
        end
        save(outfile_tex, plot, include_preamble = false)
        save(outfile_pdf, plot)
    end
end


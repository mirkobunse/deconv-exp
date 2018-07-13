"""
    smoothing(metricsfile="res/metrics/smoothing.csv")

Generate smoothing plots from experimental results.
"""
function smoothing(metricsfile::String="res/metrics/smoothing.csv"; full=false)
    initialize("smoothing", metricsfile) # initialize PGFPlots
    df = readtable(metricsfile) # read metrics
    
    # select relevant strategies and iterations
    if !full
        df = df[.&(df[:k] .<= 16, df[:k] .!= 0), :]
        valid_names = [ "no smoothing", map(o -> "polynomial (order $o)", [2, 3, 6, 12])... ]
        df = df[map(n -> in(n, valid_names), df[:name]), :]
    end
    
    # iterate metrics
    for ref in [ "oob", "ib" ], mkey in keys(METRICS)
        metric = Symbol("$(ref)_$(mkey)")
        
        # aggregate mean and quantiles (5% and 95%)
        agg = aggregate_bootstrap(df, [:name, :k], metric)
        
        # file paths
        outfile_pdf = pdfpath(metricsfile, "_" * string(metric))
        outfile_tex = texpath(metricsfile, "_" * string(metric))
        info("Plotting to $(outfile_pdf) and $(outfile_tex)")
        
        # axis with one layer per smoothing order
        plot = Axis( xlabel = "{\\sc Dsea} iteration \$k\$",
                     ylabel = METRIC_NAMES[mkey] * " to \$\\mathbf{f}\$",
                     style  = "progress axis, outside legend" )
        for (i, name) in enumerate(unique(df[:name]))
            sdf = agg[agg[:name] .== name, :]
            push!(plot, Plots.Linear( sdf[:k], sdf[:y],
                                      errorBars   = ErrorBars( yplus  = sdf[:y_plus],
                                                               yminus = sdf[:y_minus]),
                                      style       = progress_style(i),
                                      legendentry = name ))
        end
        save(outfile_tex, plot, include_preamble = false)
        save(outfile_pdf, plot)
    end
end


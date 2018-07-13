"""
    expand_reduce(metricsfile="res/metrics/expand_reduce.csv")

Generate expansion plots from experimental results.
"""
function expand_reduce(metricsfile::String="res/metrics/expand_reduce.csv"; full::Bool=false)
    initialize("expand_reduce", metricsfile) # initialize PGFPlots
    df = readtable(metricsfile) # read metrics
    
    # select relevant strategies and iterations
    if !full
        df = df[.&(df[:k] .<= 16, df[:k] .!= 0), :]
        valid_names = [ "no expansion", "expand with factor 2", "expand with factor 6" ]
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
        
        # axis with one layer per expansion factor
        plot = Axis( xlabel = "{\\sc Dsea} iteration \$k\$",
                     ylabel = METRIC_NAMES[mkey] * " to \$\\mathbf{f}\$",
                     style  = "progress axis, outside legend, ymode = log" )
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


"""
    smart_control(metricsfile="res/metrics/smart_control.csv")


"""
function smart_control(metricsfile::String="res/metrics/smart_control.csv")
    initialize("smart_control", metricsfile) # initialize PGFPlots
    df = readtable(metricsfile) # read metrics
    
    # only plot the estimates
    df = df[.|(df[:name] .== Job.SMART_CONTROL_UNIFORM,
               df[:name] .== Job.SMART_CONTROL_REDISTRICTING), :]
    
    # plot
    for ref in [ "oob", "ib" ], mkey in keys(METRICS)
        metric = Symbol("$(ref)_$(mkey)")
        
        # aggregate mean and quantiles (5% and 95%)
        agg = aggregate_bootstrap(df, [:name, :n], metric)
        
        # file paths
        outfile_pdf = pdfpath(metricsfile, "_" * string(metric))
        outfile_tex = texpath(metricsfile, "_" * string(metric))
        info("Plotting to $(outfile_pdf) and $(outfile_tex)")
        
        # axis with one layer per control strategy
        plot = Axis( xlabel = "number \$N'\$ of training examples",
                     ylabel = METRIC_NAMES[mkey] * " to \$\\mathbf{f}\$",
                     style  = join([ "progress axis",
                                     "scale = .8",
                                     "ymode = log",
                                     "xmode = log" ], ", ") )
        for (i, name) in enumerate(unique(df[:name]))
            sdf = agg[agg[:name] .== name, :]
            push!(plot, Plots.Linear( sdf[:n], sdf[:y],
                                      errorBars   = ErrorBars( yplus  = sdf[:y_plus],
                                                               yminus = sdf[:y_minus]),
                                      style       = progress_style(i),
                                      legendentry = name ))
        end
        save(outfile_tex, plot, include_preamble = false)
        save(outfile_pdf, plot)
    end
end


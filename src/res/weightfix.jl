"""
    weightfix(metricsfile="res/metrics/weightfix.csv")

Generate weightfix plots from experimental results.
"""
function weightfix(metricsfile::String="res/metrics/weightfix.csv"; full::Bool=false)
    initialize("weightfix", metricsfile) # initialize PGFPlots
    df = readtable(metricsfile) # read metrics
    
    # select interesting iterations
    df = df[.&(df[:k] .<= 12, df[:k] .!= 0), :]
    
    # iterate metrics
    for ref in [ "oob", "ib" ], mkey in keys(METRICS)
        metric = Symbol("$(ref)_$(mkey)")
        
        # aggregate mean and quantiles (5% and 95%)
        agg_full = aggregate_bootstrap(df, [:name, :dataset, :f_train, :k], metric)
        
        # iterate data sets and corresponding values of f_train
        for dataset in unique(agg_full[:dataset])
            agg_dataset = agg_full[agg_full[:dataset] .== dataset, :]
            
            for f_train in unique(agg_dataset[:f_train])
                agg = agg_dataset[agg_dataset[:f_train] .== f_train, :]
            
                # file paths
                outfile_pdf = pdfpath(metricsfile, "_" * join([dataset, f_train, string(metric)], "_"))
                outfile_tex = texpath(metricsfile, "_" * join([dataset, f_train, string(metric)], "_"))
                info("Plotting to $(outfile_pdf) and $(outfile_tex)")
                
                # axis with one layer for original DSEA and one layer for fixed weighting
                plot = Axis( xlabel = "{\\sc Dsea} iteration \$k\$",
                             ylabel = METRIC_NAMES[mkey] * " to \$\\mathbf{f}\$",
                             style  = "progress axis, outside legend" )
                for (i_name, name) in enumerate(unique(agg[:name]))
                    sdf = agg[agg[:name] .== name, :]
                    push!(plot, Plots.Linear( sdf[:k], sdf[:y],
                                              errorBars   = ErrorBars( yplus  = sdf[:y_plus],
                                                                       yminus = sdf[:y_minus]),
                                              style       = progress_style(i_name),
                                              legendentry = name ))
                end
                save(outfile_tex, plot, include_preamble = false)
                save(outfile_pdf, plot)
            end
        end
    end
end


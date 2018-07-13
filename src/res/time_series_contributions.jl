"""
    time_series_contributions(contribfile="res/time_series_contributions.csv")

Plot a simple time series of the mode of windowed deconvolution results.
"""
function time_series_contributions(contribfile::String="res/time_series_contributions.csv")
    initialize("time_series_contributions", contribfile) # initialize PGFPlots
    df = readtable(contribfile) # read metrics
    
    # file paths
    outfile_pdf = pdfpath(replace(contribfile, "res/", "res/metrics/"))
    outfile_tex = texpath(replace(contribfile, "res/", "res/metrics/"))
    info("Plotting to $(outfile_pdf) and $(outfile_tex)")
    
    # plot
    plot = Axis( Plots.Linear(df[:time], df[:mode],
                              style = joinstyles( progress_style(1),
                                                  "mark = none",
                                                  "opacity = 1" )),
                 xlabel = "index of sliding window (time)",
                 ylabel = "mode of \$\\hat{\\mathbf{f}}\$",
                 style  = "progress axis" )
    save(outfile_tex, plot, include_preamble = false)
    save(outfile_pdf, plot)
end


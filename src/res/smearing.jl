"""
    smearing(metricsfile="res/metrics/smearing.csv")

Generate smearing plots from experimental results.
"""
function smearing(metricsfile::String="res/metrics/smearing.csv")
    initialize("smearing", metricsfile) # initialize PGFPlots
    df = readtable(metricsfile) # read metrics
    
    # select interesting iterations and order results
    df = df[.&(df[:k] .<= 16, df[:k] .!= 0), :]
    order = Dict("medium smearing" => 2,
                 "weak smearing"   => 1,
                 "strong smearing" => 3)
    df[:order] = map(f -> order[f], df[:name])
    sort!(df, cols=[:order])
    
    # plot
    for ref in [ "oob", "ib" ], mkey in keys(METRICS)
        metric = Symbol("$(ref)_$(mkey)")
        
        # aggregate mean and quantiles (5% and 95%)
        agg = aggregate_bootstrap(df, [:name, :k], metric)
        
        # file paths
        outfile_pdf = pdfpath(metricsfile, "_" * string(metric))
        outfile_tex = texpath(metricsfile, "_" * string(metric))
        info("Plotting to $(outfile_pdf) and $(outfile_tex)")
        
        # axis with one layer per smearing amount
        plot = Axis( xlabel = "{\\sc Dsea} iteration \$k\$",
                     ylabel = METRIC_NAMES[mkey] * " to \$\\mathbf{f}\$",
                     style  = "progress axis" )
        for (i_name, name) in enumerate(unique(agg[:name]))
            sdf = agg[agg[:name] .== name, :]
            push!(plot, Plots.Linear( sdf[:k], sdf[:y],
                                      errorBars   = ErrorBars( yplus  = sdf[:y_plus],
                                                               yminus = sdf[:y_minus]),
                                      style       = progress_style(i_name == 1 ? 3 : i_name-1),
                                      legendentry = name ))
        end
        save(outfile_tex, plot, include_preamble = false)
        save(outfile_pdf, plot)
    end
end


_smearing_HISTOGRAM_CONFIGFILE = "conf/data/gaussian.yml" # amount of smearing
_smearing_HISTOGRAM_B          = 1                        # bootstrap iteration
_smearing_HISTOGRAM_METRIC     = :emd                     # key of Res.METRICS

"""
    smearing_histogram(spectrafile="res/spectra/smearing.csv",
                       metricsfile="res/metrics/smearing.csv")

Plot DSEA histograms of first and best iteration together with the true spectrum.
"""
function smearing_histogram(spectrafile::String="res/spectra/smearing.csv",
                            metricsfile::String="res/metrics/smearing.csv")
    
    initialize("smearing_histogram", spectrafile, metricsfile) # initialize PGFPlots
    
    # read files
    met  = readtable(metricsfile)         # metrics
    spec = unique(readtable(spectrafile)) # spectra
    
    # # find the best iteration (in-bag) in smearing and bootstrap iteration
    met = met[.&( met[:configfile] .== _smearing_HISTOGRAM_CONFIGFILE,
                  met[:b]          .== _smearing_HISTOGRAM_B,
                  met[:k]          .>= 1 ), :]
    metric = Symbol("ib_$(_smearing_HISTOGRAM_METRIC)")
    best_k = met[findmin(met[metric])[2], :k]
    info("The iteration with the best $(METRIC_NAMES[_smearing_HISTOGRAM_METRIC]) is $(best_k)")
    
    # parse deconvolution result from best_k
    best_f = eval(parse(spec[.&( spec[:configfile] .== _smearing_HISTOGRAM_CONFIGFILE,
                                 spec[:b]          .== _smearing_HISTOGRAM_B,
                                 spec[:k]          .== best_k ), :f][1]))
    f_oob  = eval(parse(spec[.&( spec[:configfile] .== _smearing_HISTOGRAM_CONFIGFILE,
                                 spec[:b]          .== _smearing_HISTOGRAM_B,
                                 spec[:name]       .== Job.TRUESPEC_OOB ), :f][1]))
    bins = binedges(Data.discretizer(Data.Gaussian(_smearing_HISTOGRAM_CONFIGFILE)))[1:end-1]
    
    # file paths
    outfile_pdf = pdfpath(metricsfile, "_histogram")
    outfile_tex = texpath(metricsfile, "_histogram")
    info("Plotting to $(outfile_pdf) and $(outfile_tex)")
    
    # plot histogram
    legend = [ "\$\\mathbf{f}\$", "\$\\hat{\\mathbf{f}}^{($(best_k))}\$"]
    style  = [ "black, thin, name path = truth, dash pattern={on 2.25pt off 0.75pt}",
               "tu01, thick, name path = estimate" ]
    plot = plot_histogram( bins,
                           (f_oob,  style[1], legend[1]),
                           (best_f, style[2], legend[2]),
                           xlabel = "value \$y \\in Y\$",
                           ylabel = "probability",
                           style  = "scale = .8, ymax = 1e0" )
    push!(plot, Plots.Command("\\addplot[draw=none, fill=tu01, fill opacity=0.33] fill between[of = truth and estimate]"))
    save(outfile_tex, plot, include_preamble = false)
    save(outfile_pdf, plot)
    
    # distance progress plot
    outfile_pdf = pdfpath(metricsfile, "_histogram_distance")
    outfile_tex = texpath(metricsfile, "_histogram_distance")
    info("Plotting to $(outfile_pdf) and $(outfile_tex)")
    plot = Axis( Plots.Linear(met[:k], met[metric], style = progress_style(1)),
                 xlabel = "{\\sc Dsea} iteration \$k\$",
                 ylabel = METRIC_NAMES[_smearing_HISTOGRAM_METRIC] * " to \$\\mathbf{f}\$",
                 style  = "progress axis, scale = .8" )
    save(outfile_tex, plot, include_preamble = false)
    save(outfile_pdf, plot)
    
end


COMPARISON_NAMES = [ "{\\sc Dsea}",
                     "\$\\text{\\sc Dsea}^+\$",
                     "{\\sc Ibu}",
                     "\$\\mathcal{RUN}\$" ]
COMPARISON_STYLES = [ joinstyles("gray!30",
                                 "every error bar/.style = {opacity=1, gray!70}"),
                      joinstyles("tu01midlight",
                                 "postaction = {pattern=crosshatch dots, pattern color = tu01!85}",
                                 "every error bar/.style = {opacity=1, tu01}"),
                      joinstyles("tu03midlight",
                                 "postaction = {pattern=north east lines, pattern color = tu03!85}",
                                 "every error bar/.style = {opacity=1, tu03}"),
                      joinstyles("tu02midlight",
                                 "postaction = {pattern=north west lines, pattern color = tu02!85}",
                                 "every error bar/.style = {opacity=1, tu02}") ]

COMPARISON_AGG_COLS = [ :name, :method, :stepsize, :f_train_i, :chi2s ]

PARAMETERS_NAMES_RUN = Dict( "\$\\mathcal{RUN}\$ (no regularization, \$J = 21\$)"  => -1,
                             "\$\\mathcal{RUN}\$ (no regularization, \$J = 84\$)"  => 2,
                             "\$\\mathcal{RUN}\$ (no regularization, \$J = 147\$)" => 3,
                             "\$\\mathcal{RUN}\$ (no regularization, \$J = 210\$)" => 4,
                             "\$\\mathcal{RUN}\$ (expand with factor 2, \$J = 84\$)"  => 6,
                             "\$\\mathcal{RUN}\$ (expand with factor 2, \$J = 147\$)" => 7,
                             "\$\\mathcal{RUN}\$ (expand with factor 2, \$J = 210\$)" => 8,
                             "\$\\mathcal{RUN}\$ (expand with factor 6, \$J = 147\$)" => 10,
                             "\$\\mathcal{RUN}\$ (expand with factor 6, \$J = 210\$)" => 11 )
PARAMETERS_NAMES_IBU = Dict( "IBU (no smoothing, \$J = 21\$)"  => 1,
                             "IBU (no smoothing, \$J = 84\$)"  => 2,
                             "IBU (no smoothing, \$J = 147\$)" => 3,
                             "IBU (no smoothing, \$J = 210\$)" => 4,
                             "IBU (smoothing order 1, \$J = 21\$)"  => 6,
                             "IBU (smoothing order 1, \$J = 84\$)"  => 7,
                             "IBU (smoothing order 1, \$J = 147\$)" => 8,
                             "IBU (smoothing order 1, \$J = 210\$)" => 9,
                             "IBU (smoothing order 2, \$J = 21\$)"  => 11,
                             "IBU (smoothing order 2, \$J = 84\$)"  => 12,
                             "IBU (smoothing order 2, \$J = 147\$)" => 13,
                             "IBU (smoothing order 2, \$J = 210\$)" => 14,
                             "IBU (smoothing order 6, \$J = 21\$)"  => 16,
                             "IBU (smoothing order 6, \$J = 84\$)"  => 17,
                             "IBU (smoothing order 6, \$J = 147\$)" => 18,
                             "IBU (smoothing order 6, \$J = 210\$)" => 19 )

K_NAMES = [ "\$\\alpha^{(k)} = 1.0\$",
            "\$\\alpha^{(k)} = 0.6\$",
            "\$\\alpha^{(k)} = 0.3\$",
            "\$\\alpha^{(k)}_{\\mathcal{RUN}}\$ with \$\\tau = 0\$",
            "\$\\alpha^{(k)}_{\\mathcal{RUN}}\$ with \$\\tau = 10^{-06}\$",
            "\$\\alpha^{(k)}_{\\mathcal{RUN}}\$ with \$\\tau = 10^{-03}\$",
            "\$\\alpha^{(k)} = k^{0.9 - 1}\$",
            "\$\\alpha^{(k)} = k^{0.6 - 1}\$",
            "\$\\alpha^{(k)} = k^{0.3 - 1}\$",
            "\$\\alpha^{(k)} = {0.9}^{k - 1}\$",
            "\$\\alpha^{(k)} = {0.6}^{k - 1}\$",
            "\$\\alpha^{(k)} = {0.3}^{k - 1}\$" ]

PDFS_FTRAIN = [ "uniform", "auxiliary" ]
PDFS_METRIC = :emd
PDFS_COLORS = Dict( "_f_dsea_constant" => "gray",
                    "_f_dsea_run"      => "tu01",
                    "_f_ibu"           => "tu03",
                    "_f_run"           => "tu02" ) # suffix -> color

"""
    all_comparisons()

Generate plots for all comparison metrics in `res/metrics`.
"""
function all_comparisons()
    metricsdir = "res/metrics" # find files in spectra directory
    files = filter(readdir(metricsdir)) do f
        startswith(f, "comparison") && endswith(f, ".csv")
    end
    
    datasets = Dict(
        "gaussian" => [ "comparison_gaussian_appropriate.csv", "comparison_gaussian_auxiliary.csv" ],
        "fact"     => [ "comparison_fact_appropriate.csv",     "comparison_fact_uniform.csv" ],
        "magic"    => [ "comparison_magic_appropriate.csv",    "comparison_magic_uniform.csv" ]
    )
    
    # select data sets for which all spectra are present
    datasets = filter((k, v) -> all(map(f -> in(f, files), v)), datasets)
    info("About to generate plots for the following metrics files:\n  - ",
         join(map(v -> join(v, " and "), values(datasets)), "\n  - "))
    
    map(values(datasets)) do fs # generate plots for each data set
        comparison(joinpath.(metricsdir, fs)...)
    end
    info("Done")
end

"""
    comparison(metricsfile...)

Generate comparison plots from experimental results.
"""
function comparison(metricsfiles::String...)
    initialize("comparison", metricsfiles...) # initialize PGFPlots
    
    # dummy path to generate pdf and tex file paths from
    outfile_csv = _comparison_common_path([metricsfiles...])
    
    # read metrics
    info("Reading $(length(metricsfiles)) metrics files..")
    df = vcat(readtable.(metricsfiles)...)
    
    # index of training density (0 or 1, given that auxiliary == uniform)
    df[:f_train_i] = zeros(Int64, size(df, 1))
    df[.|(df[:f_train] .== "uniform", df[:f_train] .== "auxiliary"), :f_train_i] = 1
    
    # select combinations from :name, :f_train, and :chi2s, which have 20 bootstrap iterations
    counts = by(df, [:name, :f_train_i, :chi2s], sdf -> DataFrame(B = size(sdf, 1))) # count iterations
    info("$(mean(counts[:B] .== 20)*100)% of the strategies have 20 bootstrap iterations.")
    counts = counts[counts[:B] .== 20, :] # select strategies
    df = join(df, counts, kind = :semi, on = [:name, :f_train_i, :chi2s]) # select from df
    
    # iterate metrics
    for mkey in keys(METRICS)
        ib_metric  = Symbol("ib_$(mkey)")
        oob_metric = Symbol("oob_$(mkey)")
        
        # aggregate in-bag metric and split into groups
        agg_ib = aggregate_bootstrap(df, COMPARISON_AGG_COLS, ib_metric)
        i_orig = agg_ib[:name]     .== STEPSIZE_ORIGINAL_NAME # original DSEA
        i_plus = agg_ib[:stepsize] .== "run"                  # DSEA+ with adaptive step size
        i_ibu  = agg_ib[:method]   .== "ibu"
        i_run  = agg_ib[:method]   .== "run"
        
        # find the best methods with respect to in-bag metric
        agg_best = vcat(map(i_group -> begin
            
            # return best of each group for each f_train distribution
            group = agg_ib[i_group, :] # group defined by its indices
            by(group, [:f_train_i]) do sdf
                sdf[findmin(sdf[:y])[2], :]
            end
            
        end, [i_orig, i_plus, i_ibu, i_run])...)
        
        # compute OOB metric of the selected strategies
        df_best = join(df, agg_best, kind = :semi, on = [:name, :chi2s, :f_train_i]) # select from df
        agg_oob = aggregate_bootstrap(df_best, COMPARISON_AGG_COLS, oob_metric)
        agg_k   = aggregate_bootstrap(df_best, COMPARISON_AGG_COLS, :k) # iteration number
        
        # extract method names from strategy names
        extract_method = n -> begin
            if contains(n, "alpha")
                COMPARISON_NAMES[2] # DSEA+
            elseif contains(n, "RUN")
                COMPARISON_NAMES[4] # RUN
            elseif contains(n, "IBU")
                COMPARISON_NAMES[3] # IBU
            elseif contains(n, "original")
                COMPARISON_NAMES[1] # DSEA
            end
        end
        agg_oob[:m] = map(extract_method, agg_oob[:name])
        agg_k[:m]   = map(extract_method, agg_k[:name])
        
        # 
        # Metrics plots
        # 
        outfile_pdf = pdfpath(outfile_csv, "_$(mkey)")
        outfile_tex = texpath(outfile_csv, "_$(mkey)")
        info("Plotting to $(outfile_pdf) and $(outfile_tex)")
        
        # axis with one layer per strategy
        barwidth = 8 # pt
        barsep   = 1
        axisstyle = join([ "default axis",
                           "outside legend",
                           "width  = .575*\\axisdefaultwidth",
                           "height = .65*\\axisdefaultheight",
                           "enlarge x limits = 0.5",
                           "ybar = $(barsep)pt",
                           "bar width = $(barwidth)pt",
                           "ymode = log",
                           "log origin = infty",
                           "xtick = {0,1}",
                           "xticklabels = {{appropriate\\strut}, {uniform\\strut}}",
                           "xticklabel style = {font=\\scriptsize}",
                           "scaled y ticks = false" ], ", ")
        if mkey == :emd
            axisstyle = joinstyles(axisstyle, "ymin = 1e-2, ymax = 1e0")
        end
        plot = Axis( xlabel = "training density",
                     ylabel = METRIC_NAMES[mkey] * " to \$\\mathbf{f}\$",
                     style  = axisstyle )
        for (i_m, m) in enumerate(COMPARISON_NAMES)
            sdf = agg_oob[agg_oob[:m] .== m, :]
            push!(plot, Plots.Linear( sdf[:f_train_i], sdf[:y],
                                      errorBars   = ErrorBars( yplus  = sdf[:y_plus],
                                                               yminus = sdf[:y_minus]),
                                      style       = joinstyles(COMPARISON_STYLES[i_m], "mark = none, draw = none"),
                                      legendentry = m ))
            # push!(plot, Plots.Command("\\node[rotate = 90, yshift = $shift, anchor = west] at ($(sdf[j,:u]), 0.011) {\\tiny\$$(sdf[j,:params])\$}"))
        end
        save(outfile_tex, plot, include_preamble = false)
        save(outfile_pdf, plot)
        
        # 
        # Iteration number plots
        # 
        outfile_pdf = pdfpath(outfile_csv, "_$(mkey)_k")
        outfile_tex = texpath(outfile_csv, "_$(mkey)_k")
        info("Plotting to $(outfile_pdf) and $(outfile_tex)")
        
        # axis with one layer per strategy
        plot = Axis( xlabel = "training density",
                     ylabel = "Iteration number \$k\$",
                     style  = joinstyles(axisstyle, "ymin = 1e0, ymax = 1e2") )
        for (i_m, m) in enumerate(COMPARISON_NAMES)
            sdf = agg_k[agg_k[:m] .== m, :]
            push!(plot, Plots.Linear( sdf[:f_train_i], sdf[:y],
                                      errorBars   = ErrorBars( yplus  = sdf[:y_plus],
                                                               yminus = sdf[:y_minus]),
                                      style       = joinstyles(COMPARISON_STYLES[i_m], "mark = none, draw = none"),
                                      legendentry = m ))
        end
        save(outfile_tex, plot, include_preamble = false)
        save(outfile_pdf, plot)
    end
end


# produce a dummy csv path from common parts of file names (used to obtain pdf and tex file paths)
function _comparison_common_path(metricsfiles::Array{String,1})
    # zip parts of file names to be compared
    z = collect(zip(map(p -> split(p, ['_', '.']), basename.(metricsfiles))...))[1:end-1]
    
    # create file name from common parts
    b = join(map(t -> t[1], z[map(t -> all(t .== t[1]), z)]), "_") # select and join parts
    return joinpath(dirname(metricsfiles[1]), b) * ".csv" # join with directory and extension
end


"""
    all_comparison_parameters(; kwargs...)

Generate plots for all comparison metrics in `res/metrics`.
"""
function all_comparison_parameters(; kwargs...)
    metricsdir = "res/metrics" # find files in spectra directory
    files = filter(readdir(metricsdir)) do f
        startswith(f, "comparison") && endswith(f, ".csv")
    end
    info("About to generate plots for the following metrics files:\n  - ", join(files, "\n  - "))
    
    map(files) do f # generate plots for each file
        comparison_parameters(joinpath(metricsdir, f); kwargs...)
    end
    info("Done")
end

"""
    comparison_parameters(metricsfile; full=false, metric=:all)

Generate plots about the meta parameter impact in the comparison experiments.
"""
function comparison_parameters(metricsfile::String; full::Bool=false, metric::Symbol=:all)
    initialize("comparison_parameters", metricsfile) # initialize PGFPlots
    df = readtable(metricsfile) # read metrics
    
    # select relevant chi2s values
    if !full
        df = df[.|(df[:chi2s] .== 1e-3, df[:chi2s] .== 1e-6), :]
    end
    
    # iterate metrics
    metrics = metric == :all ? keys(METRICS) : [ metric ]
    for ref in [ "oob", "ib" ], mkey in metrics
        ref_metric = Symbol("$(ref)_$(mkey)")
        
        # aggregate mean and quantiles (5% and 95%)
        agg_full = aggregate_bootstrap(df, [:name, :method, :stepsize, :chi2s], ref_metric)
        
        # extract parameters from strategy (used to set the x coordinate labels)
        extract_params = n -> begin
            m = match(r"\((.*?)\)", n)
            m == nothing || length(m[1]) < 3 ? n : m[1]
        end
        
        # one plot per method
        for method in unique(agg_full[:method]), stepsize in unique(agg_full[:stepsize])
            agg = agg_full[.&(agg_full[:method]   .== method,
                              agg_full[:stepsize] .== stepsize), :]
            if size(agg, 1) == 0  continue  end # several combinations do not exist
            
            # file paths
            suffix = "_$method"
            if method == "dsea"
                suffix *= "_$stepsize"
            end
            suffix *= "_$(ref_metric)"
            outfile_pdf = pdfpath(metricsfile, suffix)
            outfile_tex = texpath(metricsfile, suffix)
            info("Plotting to $(outfile_pdf) and $(outfile_tex)")
            
            # map names to x coordinates (default: indices)
            i_names = Dict(map(tup -> (tup[2], tup[1]), enumerate(unique(agg[:name]))))
            if method == "run"
                i_names = PARAMETERS_NAMES_RUN
            elseif method == "ibu"
                i_names = PARAMETERS_NAMES_IBU
            end
            agg[:i_name] = map(n -> get(i_names, n, -1), agg[:name])
            agg = agg[agg[:i_name] .!= -1, :] # remove strategies not selected by i_names
            ticks  = join(values(i_names), ",")
            labels = join(map(n -> "{$n}", map(extract_params, keys(i_names))), ",")
            
            # axis with one layer per chi2s value
            plot = Axis( xlabel = "",
                         ylabel = METRIC_NAMES[mkey] * " to \$\\mathbf{f}\$",
                         style  = joinstyles( "default axis",
                                              "scaled y ticks = false",
                                              "xtick = {$ticks}",
                                              "xticklabels = {$labels}",
                                              "xticklabel style = {font=\\tiny, rotate=90}") )
            for (i_chi2s, chi2s) in enumerate(sort(unique(agg[:chi2s])))
                sdf = agg[agg[:chi2s] .== chi2s, :]
                push!(plot, Plots.Linear( sdf[:i_name], sdf[:y],
                                          errorBars   = ErrorBars( yplus  = sdf[:y_plus],
                                                                   yminus = sdf[:y_minus]),
                                          style       = joinstyles(progress_style(i_chi2s), "only marks"),
                                          legendentry = Util.latex_e(chi2s) ))
               
            end
            save(outfile_tex, plot, include_preamble = false)
            save(outfile_pdf, plot)
        end
    end
end



"""
    all_comparison_pdfs()

Plot all densities obtained with `comparison_pdfs`.
"""
all_comparison_pdfs() =
    for fname in [ "comparison_magic_uniform",
                   "comparison_fact_uniform",
                   "comparison_gaussian_auxiliary" ]
        comparison_pdfs("res/spectra/$fname.csv", "res/metrics/$fname.csv")
    end

"""
    comparison_pdfs(spectrafile, metricsfile)

Plot the deconvolution results obtained with the best methods.
"""
function comparison_pdfs(spectrafile::String, metricsfile::String)
    initialize("comparison_pdfs", spectrafile, metricsfile) # initialize PGFPlots
    
    # obtain bins of data set
    id_dataset = convert(String, match(r"res/metrics/comparison_(.*?)_", metricsfile)[1])
    bins = binedges(Data.discretizer(Data.dataset(id_dataset, readdata = false)))[1:end-1]
    
    # read files
    info("Reading metrics and spectra of $(id_dataset) data..")
    met  = readtable(metricsfile)         # metrics
    spec = unique(readtable(spectrafile)) # spectra
    
    # split into true spectra and estimates
    spec_oob = spec[spec[:name] .== Job.TRUESPEC_OOB, :]
    spec_est = _metrics_comparison(spec) # project iterations to convergence thresholds
    
    # select configurations which have 20 bootstrap iterations
    info("Selecting spectra from best configurations..")
    counts = by(met, [:name, :chi2s], sdf -> DataFrame(B = size(sdf, 1))) # count iterations
    info("$(mean(counts[:B] .== 20)*100)% of the strategies have 20 bootstrap iterations.")
    counts = counts[counts[:B] .== 20, :] # reduce
    met = join(met, counts, kind = :semi, on = [:name, :chi2s]) # select from DataFrame
    
    # aggregate in-bag metric and split into groups
    ib_metric  = Symbol("ib_$(PDFS_METRIC)")
    met[:f_train_i] = 0.0 # presence of column assumed by aggregate_bootstrap(..)
    agg_ib = aggregate_bootstrap(met, COMPARISON_AGG_COLS, ib_metric)
    i_orig = agg_ib[:name]     .== STEPSIZE_ORIGINAL_NAME # original DSEA
    i_plus = agg_ib[:stepsize] .== "run"                  # DSEA+ with adaptive step size
    i_ibu  = agg_ib[:method]   .== "ibu"
    i_run  = agg_ib[:method]   .== "run"
    
    # find the best methods with respect to in-bag metric
    agg_best = vcat(map(i_group -> begin
        group = agg_ib[i_group, :] # group defined by its indices
        group[findmin(group[:y])[2], :]
    end, [i_orig, i_plus, i_ibu, i_run])...)
    spec_est = join(spec_est, agg_best, kind = :semi, on = [:name, :chi2s]) # select from spectra
    
    # compute the mean estimate of each strategy
    f_oob = .+(eval.(parse.(spec_oob[:f]))...) ./ size(spec_oob, 1)
    spec_est = by(spec_est, [:name, :method, :stepsize]) do sdf
        f_est = .+(eval.(parse.(sdf[:f]))...) ./ size(sdf, 1)
        f_est = f_est[f_oob .!= 0] # limit to non-zero bins
        DataFrame(f = [ f_est ])
    end
    bins  = bins[f_oob .!= 0]
    f_oob = f_oob[f_oob .!= 0]
    
    # iterate estimates
    for i in 1:size(spec_est, 1)
        f_est = spec_est[i, :f]
        
        # file paths
        method = spec_est[i, :method]
        suffix = "_f_$method"
        if method == "dsea"
            suffix *= "_$(spec_est[i, :stepsize])"
        end
        outfile_pdf = pdfpath(metricsfile, suffix)
        outfile_tex = texpath(metricsfile, suffix)
        info("Plotting to $(outfile_pdf) and $(outfile_tex)")
        
        # method name
        extract_method = n -> begin
            if contains(n, "alpha")
                COMPARISON_NAMES[2] # DSEA+
            elseif contains(n, "RUN")
                COMPARISON_NAMES[4] # RUN
            elseif contains(n, "IBU")
                COMPARISON_NAMES[3] # IBU
            elseif contains(n, "original")
                COMPARISON_NAMES[1] # DSEA
            end
        end
        method_name = extract_method(spec_est[i, :name])
        
        # plot histogram
        color  = PDFS_COLORS[suffix]
        legend = [ "\$\\mathbf{f}\$", "\$\\hat{\\mathbf{f}}_\\text{$(method_name)}\$"]
        style  = [ "black, thin, name path = truth, dash pattern={on 2.25pt off 0.75pt}",
                   "$color, thick, name path = estimate" ]
        plot = plot_histogram( bins,
                               (f_oob, style[1], legend[1]),
                               (f_est, style[2], legend[2]),
                               xlabel = "\$\\mathcal{Y}\$",
                               ylabel = "probability",
                               style  = "scale = .65" )
        push!(plot, Plots.Command("\\addplot[draw=none, fill=$color, fill opacity=0.33] fill between[of = truth and estimate]"))
        save(outfile_tex, plot, include_preamble = false)
        save(outfile_pdf, plot)
    end
    
end

"""
    all_comparison_ks(; kwargs...)

Generate plots for all comparison metrics in `res/metrics`.
"""
function all_comparison_ks(; kwargs...)
    metricsdir = "res/metrics" # find files in spectra directory
    files = filter(readdir(metricsdir)) do f
        startswith(f, "comparison") && endswith(f, ".csv")
    end
    info("About to generate plots for the following metrics files:\n  - ", join(files, "\n  - "))
    
    map(files) do f # generate plots for each file
        comparison_k(joinpath(metricsdir, f); kwargs...)
    end
    info("Done")
end

"""
    comparison_k(metricsfile...)

Generate plots on iteration numbers in DSEA/DSEA+.
"""
function comparison_k(metricsfile::String)
    initialize("comparison_k", metricsfile) # initialize PGFPlots
    df = readtable(metricsfile) # read metrics
    
    # only regard specific DSEA+ configurations
    df = df[map(n -> in(n, K_NAMES), df[:name]), :]
    df[:f_train_i] = 0.0 # presence of column assumed by aggregate_bootstrap(..)
    
    # select combinations from :name and :chi2s, which have 20 bootstrap iterations
    counts = by(df, [:name, :chi2s], sdf -> DataFrame(B = size(sdf, 1))) # count iterations
    info("$(mean(counts[:B] .== 20)*100)% of the strategies have 20 bootstrap iterations.")
    counts = counts[counts[:B] .== 20, :]
    df = join(df, counts, kind = :semi, on = [:name, :chi2s]) # select from df

    # aggregate iteration numbers
    agg_full = aggregate_bootstrap(df, COMPARISON_AGG_COLS, :k)
    
    # iterate step size strategies
    for (i_stepsize, stepsize) in enumerate(unique(agg_full[:stepsize]))
        agg = agg_full[agg_full[:stepsize] .== stepsize, :]
        
        # file paths
        outfile_pdf = pdfpath(metricsfile, "_$(stepsize)_k")
        outfile_tex = texpath(metricsfile, "_$(stepsize)_k")
        info("Plotting to $(outfile_pdf) and $(outfile_tex)")
        
        # axis with one layer per stepsize strategy
        plot = Axis( xlabel = "convergence threshold \\epsilon",
                     ylabel = "iteration number \$k\$",
                     style  = joinstyles("progress axis",
                                         "outside legend",
                                         "xmode = log",
                                         "ymode = log",
                                         "ymin = 1e0",
                                         "ymax = 1e2",
                                         "x dir = reverse") )
        for (i_name, name) in collect(enumerate(unique(agg[:name])))[end:-1:1]
            sdf = agg[agg[:name] .== name, :]
            push!(plot, Plots.Linear( sdf[:chi2s], sdf[:y],
                                      errorBars   = ErrorBars( yplus  = sdf[:y_plus],
                                                               yminus = sdf[:y_minus]),
                                      style       = progress_style(i_name),
                                      legendentry = name ))
        end
        save(outfile_tex, plot, include_preamble = false)
        save(outfile_pdf, plot)
    end
end



STEPSIZE_ORIGINAL_NAME = "original {\\sc Dsea}"      # name of original DSEA strategy
STEPSIZE_OPTIMAL_NAME  = "optimal \$\\alpha^{(k)}\$" # name of optimal step size strategy

"""
    all_stepsizes()

Generate plots for all stepsize metrics in `res/metrics`.
"""
function all_stepsizes()
    metricsdir = "res/metrics" # find files in spectra directory
    files = filter(readdir(metricsdir)) do f
        startswith(f, "stepsize") && endswith(f, ".csv")
    end
    info("About to generate plots for the following metrics files:\n  - ", join(files, "\n  - "))
    
    map(files) do f # generate plots for each file
        stepsize(joinpath(metricsdir, f))
    end
    info("Done")
end

"""
    stepsize(metricsfile; full=false)

Generate stepsize plots from experimental results.
"""
function stepsize(metricsfile::String; full=false)
    initialize("stepsize", metricsfile) # initialize PGFPlots
    df = readtable(metricsfile) # read metrics
    
    # select relevant strategies and iterations
    if !full
        df = df[.&(df[:k] .<= 16, df[:k] .!= 0), :]
        valid_names = vcat(
            # valid constant names
            filter(n -> contains(n, "0.3") || contains(n, "0.6"),
                   unique(df[df[:method] .== "constant", :name])),
            
            # valid decay_mul names
            filter(n -> contains(n, "0.0") || contains(n, "0.5"),
                   unique(df[df[:method] .== "decay_mul", :name])),
            
            # valid decay_exp names
            filter(n -> contains(n, "0.3") || contains(n, "0.6"),
                   unique(df[df[:method] .== "decay_exp", :name])),
            
            # valid run names
            filter(n -> contains(n, "{-06}") || contains(n, "\\tau = 0"),
                   unique(df[df[:method] .== "run", :name])),
            
            # strategies that are always plotted
            STEPSIZE_ORIGINAL_NAME,
            STEPSIZE_OPTIMAL_NAME
        )
        df = df[map(n -> in(n, valid_names), df[:name]), :]
    end
    
    # iterate metrics
    for ref in [ "oob", "ib" ], mkey in keys(METRICS)
        metric = Symbol("$(ref)_$(mkey)")
        
        # aggregate mean and quantiles (5% and 95%)
        agg_full = aggregate_bootstrap(df, [:name, :method, :k], metric)
        
        # isolate methods that are always plotted
        agg_orig = agg_full[agg_full[:name] .== STEPSIZE_ORIGINAL_NAME, :]
        agg_opt  = agg_full[agg_full[:name] .== STEPSIZE_OPTIMAL_NAME,  :]
        agg_full = agg_full[.&(agg_full[:name] .!= STEPSIZE_ORIGINAL_NAME,
                               agg_full[:name] .!= STEPSIZE_OPTIMAL_NAME), :]
        
        # one plot per method
        for method in unique(agg_full[:method])
            agg = vcat(agg_orig, agg_opt, agg_full[agg_full[:method] .== method, :])
            
            # file paths
            outfile_pdf = pdfpath(metricsfile, "_" * join([method, string(metric)], "_"))
            outfile_tex = texpath(metricsfile, "_" * join([method, string(metric)], "_"))
            info("Plotting to $(outfile_pdf) and $(outfile_tex)")
            
            # axis with one layer per stepsize strategy
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

# 
# TODO Repair the following MDS plots
# 

STANDARD_NAME = "original"
OPTIMAL_NAME  = "optimal"

"""
    stepsize_mds(spectrafile)

Multidimensional Scaling plots from stepsize results.
"""
function stepsize_mds(spectrafile::String)
    initialize("stepsize_mds", spectrafile)
    
    # read data and rename strategies
    df = unique(readtable(spectrafile))
    df[df[:name] .== "\$\\alpha = 1.0\$", :name] = STANDARD_NAME
    df[df[:name] .== "optimal", :name]           = OPTIMAL_NAME
    df[df[:name] .== "\$\\alpha_\\mathcal{RUN} \\; (\\tau = 0e+00)\$", :name] = "\$\\alpha_\\mathcal{RUN} \\; (\\tau = 0)\$"
    
    _stepsize_mds(df, [
        STANDARD_NAME,
        OPTIMAL_NAME,
        "\$\\alpha_\\mathcal{RUN} \\; (\\tau = 0)\$"
    ], spectrafile, "_mds")
end

"""
    stepsize_smoothing_mds(spectrafile, smoothingspectrafile)

Multidimensional Scaling plots comparing stepsize and smoothing strategies.
"""
function stepsize_smoothing_mds(spectrafile::String, smoothingspectrafile::String)
    initialize("stepsize_smoothing_mds", spectrafile, smoothingspectrafile)
    
    # read data
    df = unique(vcat(readtable(spectrafile), readtable(smoothingspectrafile)))
    
    # define selection and convert to lowercase
    strategies = [
        "No Smoothing",
        "Polynomial (order 2)",
        "Polynomial (order 6)",
        "Polynomial (order 12)" ]
    indices = [ in(s, strategies) for s in df[:name] ]
    df[indices, :name] = lowercase.(df[indices, :name])
    strategies = lowercase.(strategies)
    
    # plot
    _stepsize_mds(df, strategies, spectrafile, "_smoothing_mds")
    
end

function _stepsize_mds(df::DataFrame, strategies::Array{String, 1}, spectrafile::String, with::String)
    
    # perform pre-selection
    df[:f] = map(s -> eval(parse(s)), df[:f]) # parse array
    df = df[(df[:k] .<= 16) & (df[:k] .!= 0), :] # select iterations
    df = df[[ in(s, vcat(strategies, ["TRUE_SPECTRUM_IN_BAG", "TRUE_SPECTRUM_OOB"])) for s in df[:name] ], :]
    
    for b in 1:3 # iterate over first three bootstraps
        
        bdf    = df[df[:b] .== b, :] # select bootstrap
        itruth = findfirst(bdf[:name] .== "TRUE_SPECTRUM_IN_BAG") # index of true spectrum
        ioob   = findfirst(bdf[:name] .== "TRUE_SPECTRUM_OOB")
        bdf    = vcat(bdf[[itruth, ioob], :],
                      bdf[[ !startswith(s, "TRUE_SPECTRUM") for s in bdf[:name] ], :])
        
        # matrix of spectra and distances between them
        spectra   = reshape(vcat(bdf[:f]...), (length(bdf[1,:f]), size(bdf,1)))
        distances = Util.pairwise_mdpa(spectra)
        
        # scatter coordinates from multi-dimensional scaling
        mds = classical_mds(distances, 2)
        bdf[:x] = mds[1, :] - mds[1, itruth] # center around truth
        bdf[:y] = mds[2, :] - mds[2, itruth]
        
        # remove true spectrum (0,0) from DataFrame
        oob = bdf[ioob, :]
        bdf = bdf[setdiff(1:size(bdf, 1), [itruth, ioob]), :]
        ifirst = findfirst(bdf[:k] .== 1) # position of first iteration
        
        # truth in lower left corner
        xrange = extrema(bdf[:x])
        yrange = extrema(bdf[:y])
        if abs(bdf[ifirst,:x] - xrange[1]) > abs(bdf[ifirst,:x] - xrange[2])
            bdf[:x] = - bdf[:x]
        end
        if abs(yrange[1]) > abs(yrange[2])
            bdf[:y] = - bdf[:y]
        end
        
        # plot
        pdfpath = _pdfpath(spectrafile, with * "_$b")
        texpath = _texpath(spectrafile, with * "_$b")
        info("Plotting to $pdfpath and $texpath")
        
        plot = Axis(Plots.Linear([0], [0], legendentry = "in-bag truth",
                                 style = "black, only marks, mark = x, mark options = {mark size = 4}"),
                    xlabel="", ylabel="",
                    style="legend pos = north east, axis equal, ticks = none, scale = .525,
                    axis line style={anthracite},
                    legend style = {draw = none, font=\\small, at = {(1.05,.5)}, anchor = west},
                    legend cell align={left}, axis background/.style = {fill = white}")
        push!(plot, Plots.Linear(oob[:x], oob[:y], legendentry = "out-of-bag truth",
                                 style = "black, only marks, mark = otimes"))
        for (i, strategy) in enumerate(strategies) # one layer per strategy
            sdf = bdf[bdf[:name] .== strategy, :]
            push!(plot, Plots.Linear(sdf[:x], sdf[:y],
                                     style=_style(i),
                                     legendentry=strategy))
        end
        push!(plot, Plots.Command("\\node[anthracite, above right = 12pt, xshift = -5pt, inner sep = 3pt] (k) at ($(bdf[ifirst,:x]), $(bdf[ifirst,:y])) {\\footnotesize\$k = 1\$}"))
        push!(plot, Plots.Command("\\draw[anthracite, densely dashed, very thin] (k) -- ($(bdf[ifirst,:x]), $(bdf[ifirst,:y]))"))
        if b < 3 # omit legend in each but the last plot
            push!(plot, Plots.Command("\\legend{}"))
        end
        save(texpath, plot, include_preamble = false)
        save(pdfpath, plot)
    
    end
end


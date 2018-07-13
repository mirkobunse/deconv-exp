"""
    data_iact()

Generate plots and tables about IACT data sets.
"""
function data_iact()
    initialize("data_iact") # initialize PGFPlots
    
    # the data sets to consider
    ids    = [ "fact", "magic" ]
    names  = uppercase.(ids) # legend entries
    styles = [ "semithick, black, name path = h1, densely dashed",
               "semithick, black, name path = h2" ] # plot styles
    
    # style of the histogram axis
    axisstyle = joinstyles( "xmode = log",
                            "outside legend",
                            "scale = .6",
                            "ylabel style = {rotate=-90, at={(axis description cs:-.25,.5)}}" )
    
    # generate table rows and histogram plots
    rows  = String[]
    plots = Axis[]
    f_min = Inf
    for (id, name, style) in zip(ids, names, styles)
        
        # read data
        info("Reading $id data..")
        dataset  = Data.dataset(id, nobs = -1) # full data set
        y_full_c = Data.y_data(dataset) # continuous values
        
        # discretizer from y_min to y_max
        y_min, y_max = extrema(y_full_c)
        discr = LinearDiscretizer(linspace(y_min, y_max, length(Data.bins(dataset))+1))
        bins  = 10 .^ binedges(discr)[1:end-1] # lower bin boundaries in log10 space
        
        # table row
        nobs = length(y_full_c)
        push!(rows, @sprintf " %9s  &  %9s \$\\gamma\$ particles  &  \$10^{%.1f}\$ GeV  &  \$10^{%.1f}\$ GeV  \\\\" name _fcomma(nobs) y_min y_max)
        
        # histogram
        f = fit_pdf(encode(discr, y_full_c), Data.bins(discr))
        f_min = min(f_min, minimum(f))
        push!(plots, plot_histogram( bins, (f,  style, name),
                                     xlabel = "\$y\$ [GeV]",
                                     ylabel = "\$\\mathbf{f}\$",
                                     style  = axisstyle,
                                     ymin   = 1e-12 ))
    end
    
    # write LaTeX table to file
    tabular = replace("""
                      % 
                      % This file was generated with Res.data_iact()
                      % 
                      % git commit = $(Git.commithash())
                      % git origin = $(Git.remoteurl())
                      % uncommited changes = $(Git.haschanges())
                      % 
                      \\begin{tabular}{rccc}
                        \\toprule
                        telescope & number of examples & \$\\min(E)\$ & \$\\max(E)\$ \\\\
                        \\midrule
                      ROWS
                        \\bottomrule
                      \\end{tabular}
                      """, "ROWS", join(rows, "\n"))
    outfile_tab = "res/table-tex/data_iact.tex"
    info("Writing to $outfile_tab..")
    fio = open(outfile_tab, "w")
    println(fio, tabular)
    close(fio)
    println(tabular) # also print to STDOUT
    
    # file paths of plots
    outfile_pdf = "res/pdf/data_iact.pdf"
    outfile_tex = "res/tex/data_iact.tex"
    info("Plotting to $outfile_pdf and $outfile_tex")
    
    # combine plots
    plot = plots[1]
    push!(plot, plots[2].plots...)
    plot.xmin = min(plot.xmin, plots[2].xmin)
    plot.xmax = max(plot.xmax, plots[2].xmax)
    plot.ymin = 10 ^ (log10(f_min) - .25)
    for i in 1:length(plots)
        plot.plots[i].legendentry = names[i]
        xmin, xmax = extrema(plot.plots[i].data[1,:])
        push!(plot, Plots.Command("\\addplot[draw = none, domain = $xmin:$xmax, name path = b$i] {$(plot.ymin)}"))
        pattern = i == 1 ? "north east lines" : "north west lines"
        push!(plot, Plots.Command("\\addplot[draw = none, pattern color=tu0$(i)midlight, pattern = $pattern] fill between[of = b$i and h$i]"))
    end
    save(outfile_tex, plot, include_preamble = false)
    save(outfile_pdf, plot)
    
end

"""
    data_toy()

Generate plots about observables in the gaussian data set.
"""
function data_toy()
    initialize("data_toy", font="\\footnotesize") # initialize PGFPlots
    
    # read data
    info("Reading data..")
    dataset = Data.Gaussian()
    X_data  = Data.X_data(dataset)
    X_train = Data.X_train(dataset)
    
    # plot distribution of first three observables
    for i in 1:3
        x_data_c  = X_data[:,  i] # continuous observable
        x_train_c = X_train[:, i]
        
        # discretizer for the observable
        x_min, x_max = extrema(vcat(x_data_c, x_train_c))
        discr = LinearDiscretizer(linspace(x_min, x_max, length(Data.bins(dataset))+1))
        bins  = binedges(discr)[1:end-1] # lower bin edges
        
        # obtain observable distributions
        g_data  = fit_pdf(encode(discr, x_data_c),  1:length(bins))
        g_train = fit_pdf(encode(discr, x_train_c), 1:length(bins))
        
        # output file paths
        outfile_pdf = "res/pdf/data_toy_$i.pdf"
        outfile_tex = "res/tex/data_toy_$i.tex"
        info("Plotting to $outfile_pdf and $outfile_tex")
        
        # plot histograms
        plot  = plot_histogram( bins,
                                (g_train, "semithick, black, densely dashed", "training data"),
                                (g_data,  "semithick, black", "observed data"),
                                xlabel = "\$x \\in X_$i\$",
                                ylabel = "probability",
                                style  = "scale = .4, outside legend, ymax = 2e0",
                                ymin   = 1e-6 )
        if i < 3
            push!(plot, Plots.Command("\\legend{}")) # omit legend in all plots but the last
        end
        if i > 1
            plot.ylabel = "" # omit y label in all plots but the first
            plot.style  = plot.style * ", yticklabels=\\empty"
        end
        save(outfile_tex, plot, include_preamble = false)
        save(outfile_pdf, plot)
        
    end
    
end


EXP_N_OBS = 18394
function data_quality_f()
    initialize("data_quality_f")          # initialize PGFPlots
    ids = [ "gaussian", "magic", "fact" ] # data sets
    srand(1337)                           # reproducible splits
    
    rows = String[]
    plot = Axis( xlabel = "data set size \$N\$",
                 ylabel = "{\\sc Emd} to \$\\mathbf{f}\$",
                 style  = "progress axis, outside legend, scale = .5, ymode = log, xmode = log" )
    
    # add one layer per data set
    for (i, id) in enumerate(ids)
        
        # read data
        info("Reading $id data..")
        dataset = Data.dataset(id, nobs = -1)
        discr   = Data.discretizer(dataset)
        bins    = Data.bins(discr)
        y_data  = encode(discr, Data.y_data(dataset))
        
        # obtain the true histogram
        f_true = fit_pdf(y_data, bins)
        
        # experiment
        info("Increasing N in multiple permutations...")
        perms = [ randperm(length(y_data)) for _ in 1:200 ] # 200 random permutations
        df    = DataFrame(N = Int64[], p = Int64[], emd = Float64[]) # result frame
        for N in vcat(convert.(Int64, round.(logspace(2, log10(length(y_data)), 10)))[1:end-1], EXP_N_OBS)
            for (p, perm) in enumerate(perms)
                f_N = fit_pdf(y_data[perm[1:N]], bins)
                emd = Util.mdpa(f_N, f_true)
                push!(df, [ N, p, emd ])
            end
            println("N = $N")
        end
        
        # aggregate results
        agg = by(df, :N) do sdf
            mean   = Base.mean(sdf[:emd])
            lo, hi = quantile(sdf[:emd], [.05, .95])
            DataFrame(mean  = mean,
                      plus  = hi - mean,
                      minus = mean - lo,
                      std   = std(sdf[:emd], mean = mean))
        end
        
        # fill table
        name = id == "gaussian" ? "Toy Data" : uppercase(id)
        N = findfirst(agg[:N] .== EXP_N_OBS)
        emd_mean, emd_std = Util.latex_e.((agg[N, :mean], agg[N, :std]), 2, dollars = false)
        push!(rows, @sprintf " %9s  &  \$%s +- %s\$ \\\\" name emd_mean emd_std)
        
        # plot
        agg = agg[agg[:N] .!= EXP_N_OBS, :]
        push!(plot, Plots.Linear( agg[:N], agg[:mean],
                                  errorBars   = ErrorBars(yplus = agg[:plus], yminus = agg[:minus]),
                                  style       = progress_style(i),
                                  legendentry = name ))
    end
    push!(plot, Plots.Command("\\draw ({axis cs:$EXP_N_OBS,0}|-{rel axis cs:0,0}) -- ({axis cs:$EXP_N_OBS,0}|-{rel axis cs:0,1})"))
    
    # write LaTeX table to file
    tabular = replace("""
                      % 
                      % This file was generated with Res.data_quality_f(\"$(join(ids, "\", \""))\")
                      % 
                      % git commit = $(Git.commithash())
                      % git origin = $(Git.remoteurl())
                      % uncommited changes = $(Git.haschanges())
                      % 
                      \\begin{tabular}{rl}
                        \\toprule
                        data set  &  {\\sc Emd} to \$\\mathbf{f}\$  \\\\
                        \\midrule
                      ROWS
                        \\bottomrule
                      \\end{tabular}
                      """, "ROWS", join(rows, "\n"))
    outpath = "res/table-tex/data_quality_f.tex"
    info("Writing to $outpath..")
    fio = open(outpath, "w")
    println(fio, tabular)
    close(fio)
    println(tabular) # also print to STDOUT
    
    # save plot
    outfile_pdf = "res/pdf/data_quality_f.pdf"
    outfile_tex = "res/tex/data_quality_f.tex"
    info("Plotting to $outfile_pdf and $outfile_tex")
    save(outfile_tex, plot, include_preamble = false)
    save(outfile_pdf, plot)
    
end
    

# utility formating a big number string for LaTeX tables
function _fcomma(n::Int)
    sgn = sign(n)
    str = string(abs(n))[end:-1:1]
    if length(str) > 3
        q = collect(3:3:(length(str)-1))[end:-1:1]
        while length(q) > 0
            i = pop!(q)
            str = replace(str, str[1:i], str[1:i] * ",", 1)
            q = q .+ 1
        end
    end
    if sgn < 0
        str = str * "-"
    end
    return str[end:-1:1]
end


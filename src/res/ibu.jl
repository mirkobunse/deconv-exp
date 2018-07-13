"""
    ibu_progress()

Plot the first and fifth iteration of the Iterative Bayesian Unfolding
"""
function ibu_progress() # function must not have the same name as CherenkovDeconvolution.ibu
    initialize("ibu_progress") # initialize PGFPlots
    
    # read data
    dataset  = Data.Gaussian()
    discr    = Data.discretizer(dataset)
    bins     = Data.bins(discr)
    X_data   = Data.X_data(dataset)
    X_train  = Data.X_train(dataset)
    y_data   = encode(discr, Data.y_data(dataset))
    y_train  = encode(discr, Data.y_train(dataset))
    
    # discretize the feature space with a decision tree from CherenkovDeconvolution.Sklearn
    discr_x = TreeDiscretizer(X_train, y_train, 210) # 210 observable bins
    x_data  = encode(discr_x, X_data)
    x_train = encode(discr_x, X_train)
    
    # Iterative Bayesian Unfolding
    fs = Array{Float64,1}[] # array of arrays
    ibu(x_data, x_train, y_train, bins, K = 5, inspect = (f, k, chi2s) -> push!(fs, f))
    
    # plot
    f_true = ( fit_pdf(y_data, bins),
               "black, thin, dash pattern={on 2.25pt off 0.75pt}, name path = truth",
               "\$\\mathbf{f}\$" ) # tuple of the truth
    cols = Dict( 1 => "tu02", 5 => "tu01" ) # color map
    for (i, col) in cols
        outfile_pdf = "res/pdf/ibu_progress_$i.pdf"
        outfile_tex = "res/tex/ibu_progress_$i.tex"
        info("Plotting to $outfile_pdf and $outfile_tex")
        
        plot  = plot_histogram( bins, f_true,
                                ( fs[i+1], "$col, thick, name path = estimate",
                                  "\$\\hat{\\mathbf{f}}^{($i)}\$" ),
                                xlabel = "value \$y \\in Y\$",
                                ylabel = "probability",
                                style  = "scale = .65, ymax = 1e0",
                                ymin   = 5e-5 )
        push!(plot, Plots.Command("\\addplot[draw = none, fill = $col, fill opacity=0.333] fill between[of = truth and estimate]"))
        save(outfile_tex, plot, include_preamble = false)
        save(outfile_pdf, plot)
    end
end


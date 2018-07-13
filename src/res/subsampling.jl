"""
    subsampling(ids = ["magic", "fact", "gaussian"])

Generate table about uniform subsampling.
"""
function subsampling(ids::AbstractArray{String,1} = ["magic", "fact", "gaussian"])
    srand(1337) # reproducible splits
    
    # generate table rows
    rows = String[]
    for (j, id) in enumerate(ids)
        name = id == "gaussian" ? "Toy Data" : uppercase(id)
        f_train_arr = [ "appropriate", id == "gaussian" ? "auxiliary" : "uniform" ]
        
        # read data
        dataset = Data.dataset(id, nobs = -1, nobs_train = -1) # full data set
        discr   = Data.discretizer(dataset)
        bins    = Data.bins(discr)
        X_full  = Data.X_data(dataset)
        y_full  = encode(discr, Data.y_data(dataset))
        
        X_auxtrain, y_auxtrain = try
            (Data.X_train(dataset), encode(discr, Data.y_train(dataset)))
        catch (nothing, nothing) end
        
        # split with different f_train configurations, obtaining the table rows
        for f_train in f_train_arr
            _, _, _, y_train = Util.shuffle_split_subsample( X_full, y_full,
                                                             X_auxtrain, y_auxtrain,
                                                             f_train,
                                                             nobs_train = -1 )
            nobs = length(y_train) # only number of training examples is relevant
            push!(rows, @sprintf " %9s  &  %16s  &  %13s ex.  \\\\" name f_train _fcomma(nobs))
        end
        
        # group data sets
        if j < length(ids)
            push!(rows, "  \\midrule")
        end
    end
    
    # write LaTeX table to file
    tabular = replace("""
                      % 
                      % This file was generated with Res.subsampling(\"$(join(ids, "\", \""))\")
                      % 
                      % git commit = $(Git.commithash())
                      % git origin = $(Git.remoteurl())
                      % uncommited changes = $(Git.haschanges())
                      % 
                      \\begin{tabular}{rll}
                        \\toprule
                        data set  &  training density  &  training pool size  \\\\
                        \\midrule
                      ROWS
                        \\bottomrule
                      \\end{tabular}
                      """, "ROWS", join(rows, "\n"))
    outpath = "res/table-tex/subsampling.tex"
    info("Writing to $outpath..")
    fio = open(outpath, "w")
    println(fio, tabular)
    close(fio)
    println(tabular) # also print to STDOUT
    
end


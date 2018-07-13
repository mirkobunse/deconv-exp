"""
    time_series_contributions(configfile="conf/job/time_series_contributions.yml")


"""
function time_series_contributions(configfile::String="conf/job/time_series_contributions.yml")
    
    # read configuration
    c = load_file(configfile)
    srand(c["seed"])
    dataset_id = c["dataset"]["id"]
    skconfig   = joinpath("conf/skl", c["dataset"]["skconfig"] * ".yml")
    window_width = c["window"]["width"]
    window_step  = c["window"]["step"]
    outfile = c["outfile"]
    
    # read data
    dataset = Data.dataset(dataset_id, nobs = -1)
    discr   = Data.discretizer(dataset)
    bins    = Data.bins(discr)
    X_full  = Data.X_data(dataset)
    y_full  = encode(discr, Data.y_data(dataset))
    
    # split into training and observed data sets (do not subsample!)
    X_data, y_data, X_train, y_train = Util.shuffle_split_subsample( X_full, y_full,
                                                                     nothing, nothing,
                                                                     "appropriate",
                                                                     nobs_train = -1 )
    
    # order observed data in low and high energy groups
    i_low  = y_data .< bins[convert(Int64, round(length(bins) / 2))] # smaller than center bin
    X_data = vcat(X_data[i_low, :], X_data[.!(i_low), :])
    y_data = vcat(y_data[i_low],    y_data[.!(i_low)])
    
    # deconvolve
    info("Deconvolving $dataset_id data..")
    tp = train_and_predict_proba(Util.classifier_from_config(skconfig))
    _, contributions = dsea(X_data, X_train, y_train, tp, bins, return_contributions = true)
    
    # compute the estimate in each window
    info("Windowing..")
    spectra = map(1:window_step:(size(contributions, 1) - window_width)) do i
        CherenkovDeconvolution._dsea_reconstruct(contributions[i:(i+window_width-1), :])
    end
    
    # create a DataFrame from artificial time stamps and the mode of each spectrum
    df = DataFrame( time = 1:length(spectra),
                    mode = map(f -> bincenters(discr)[findmax(f)[2]], spectra) )
    
    # output
    writetable(outfile, df)
    info("Results written to $outfile")
    df
    
end


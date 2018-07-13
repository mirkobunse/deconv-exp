module Res


using DataFrames, Query, Distances, Discretizers, MultivariateStats, Colors

using ComfyCommons: info, warn, Git, Pgfplots
using CherenkovDeconvolution: Util, Sklearn, ibu

import Util, Data, Job


"""
    initialize(job, files...[; font = "\\small"])

Initialize the PGFPlots backend for the function Job.`job` called for the `files`.
"""
function initialize(job::String, files::String...; font="\\small")
    
    # warn about changes
    commit     = Git.commithash()
    remote     = Git.remoteurl()
    haschanges = Git.haschanges("src", files...)
    if haschanges
        warn("Uncommited changes may affect plots")
    end
    
    # set preamble for PDF output
    preamble = join(readlines("src/res/tex/preamble.tex", chomp=false)) # basic preamble
    replacements = Dict("GitCommit={}" =>            "GitCommit={$commit}",
                        "GitOrigin={}" =>            "GitOrigin={$remote}",
                        "GitUncommitedChanges={}" => "GitUncommitedChanges={$haschanges}")
    for (k, v) in replacements
        preamble = replace(preamble, k, v)
    end
    resetPGFPlotsPreamble()
    pushPGFPlotsPreamble(preamble)
    
    # set file template
    Pgfplots.setpgftemplate("""
    % 
    % This file was generated with Res.$job(\"$(join(files, "\""))\")
    % 
    % git commit = $commit
    % git origin = $remote
    % uncommited changes = $haschanges
    % 
    \\tikzstyle{every node} = [font=$font]

    $(TEMPLATE_PLOT_PLACEHOLDER)
    """)
    
end


METRICS = Dict( :kl    => Util.kl,
                :h     => Distances.hellinger,
                :chi2p => Util.chi2p,
                :chi2s => Util.chi2s,
                :emd   => Util.mdpa )

METRIC_NAMES = Dict( :kl    => "Kullback-Leibler distance",
                     :h     => "Hellinger distance",
                     :chi2p => "\$\\mathcal{X}^2_P\$ distance",
                     :chi2s => "\$\\mathcal{X}^2_{Sym}\$ distance",
                     :emd   => "{\\sc Emd}" )

"""
    all_metrics()

Generate metrics files for all results in `res/spectra`.
"""
function all_metrics()
    spectradir = "res/spectra" # find files in spectra directory
    files = filter(readdir(spectradir)) do f
        !startswith(f, "TEST-") && endswith(f, ".csv")
    end
    info("About to compute metrics for the following spectra files:\n  - ", join(files, "\n  - "))
    
    map(files) do f # compute metrics for each file
        metrics(joinpath(spectradir, f))
    end
end

"""
    metrics(spectrafile)

Evaluate the deconvolution results stored in the `spectrafile`, writing to a corresponding
metrics file.
"""
function metrics(spectrafile::String)
    jobname = _job_name(spectrafile)
    
    # read file
    info("Reading $spectrafile..")
    df = readtable(spectrafile)
    
    # split into reference spectra and deconvolution results
    refs = Dict( "ib"  => df[df[:name] .== Job.TRUESPEC_IB,  :],
                 "oob" => df[df[:name] .== Job.TRUESPEC_OOB, :] )
    for (id, ref) in refs
        ref = unique(ref)
        ref[:f] = [ normalizepdf(eval.(parse.(ref[:f]))...)... ]
        refs[id] = ref
    end
    df = df[.&(df[:name] .!= Job.TRUESPEC_IB,
               df[:name] .!= Job.TRUESPEC_OOB), :]
    
    # special pre-processing of results
    if jobname == "comparison" # other experiments do not need this
        info("Selecting relevant results..")
        df = _metrics_comparison(df)
    end
    df = unique(df)
    
    # parse and normalize arrays
    info("Parsing $(size(df, 1)) spectra..")
    df[:f] = [ normalizepdf(eval.(parse.(df[:f]))...)... ] # may take some time
    
    # columns to join on
    jcols = get(JCOLS, jobname, [:b]) # join on :b, by default
    info("Joining results on $jcols")
    
    # evaluate metrics with respect to both reference spectra
    for (id, ref) in refs
        
        # join results with reference
        rename!(ref, :f, :f_ref)
        df = join(df, ref[:, vcat(jcols, :f_ref)], on = jcols) # join on jcols
        
        # evaluate metrics
        for (mkey, mfun) in METRICS
            mcol     = Symbol("$(id)_$(string(mkey))")
            df[mcol] = map(kv -> mfun(kv...), zip(df[:f], df[:f_ref]))
        end
        
        # remove reference spectra
        df = df[:, setdiff(names(df), [:f_ref])]
        
    end
    
    # remove result spectra
    df = df[:, setdiff(names(df), [:f])]
    
    # store metrics in file
    outfile = replace(spectrafile, "/spectra/", "/metrics/")
    info("Writing metrics to $outfile..")
    writetable(outfile, df)
    
end

_metrics_comparison(df::DataFrame) =
    vcat(map( c -> begin
    
        info("Let chi2s = $c..")
        @from i in df begin
            @where abs(i.chi2s) <= c # only consider methods that converged
            @group i by (i.name, i.f_train, i.b) into g
            @let i_min = findmin(g..k)[2] # first k with chi2s <= c
            @select { name           = (g..name)[i_min],
                      method         = (g..method)[i_min],
                      fixweighting   = (g..fixweighting)[i_min],
                      stepsize       = (g..stepsize)[i_min],
                      smoothing      = (g..smoothing)[i_min],
                      discretization = (g..discretization)[i_min],
                      J              = (g..J)[i_min],
                      f_train        = (g..f_train)[i_min],
                      b              = (g..b)[i_min],
                      k              = (g..k)[i_min],
                      chi2s          = c,
                      f              = (g..f)[i_min] }
            @collect DataFrame
        end
        
    end, 10.0 .^ -(1:9))...) # 1e-1, 1e-2, ..., 1e-9

# columns to join reference spectra on (default is [ :b ])
JCOLS = Dict( "smearing"   => [ :b, :configfile ],
              "weightfix"  => [ :b, :dataset, :f_train ],
              "comparison" => [ :b, :f_train ] )

# name of job that produced the result file (used to obtain the keys for JCOLS)
_job_name(spectrafile::String) = split(basename(spectrafile), [ '_', '.' ])[1]


"""
    aggregate_bootstrap(df, split_by, metric[; quantiles = (.05, .95)])

Aggregate the `metric` from all bootstrap samples in the DataFrame `df`, which is split into
groups with the columns in `split_by`.

Returns the mean value and the lower and upper quantiles in another DataFrame with the
columns `y`, `y_plus`, and `y_minus`.
"""
aggregate_bootstrap(df::DataFrame, split_by::AbstractArray{Symbol,1}, metric::Symbol;
                    quantiles::Tuple{Float64,Float64}=(.05, .95)) =
    by(df, split_by) do sdf                        # split into aggregation groups
        m = sdf[isfinite.(sdf[metric]), metric]    # metric array without NaNs and Infs
        m_mean     = mean(m)                       # mean value
        m_lo, m_hi = quantile(m, [ quantiles... ]) # lower and upper quantiles
        
        # return value is the aggregation of this group
        DataFrame( y       = m_mean,
                   y_plus  = m_hi - m_mean,
                   y_minus = m_mean - m_lo )
    end


"""
    pdfpath(path, suffix="")

Obtain a file path to store a `.pdf` file in, which corresponds to the input file `path`.
"""
pdfpath(path::String, suffix::String="") = _outfilepath(path, "pdf", suffix)

"""
    texpath(path, suffix="")

Obtain a file path to store a `.tex` file in, which corresponds to the input file `path`.
"""
texpath(path::String, suffix::String="") = _outfilepath(path, "tex", suffix)

_outfilepath(path::String, ext::String, suffix::String) =
    replace(replace(path, "/metrics/", "/$ext/"), r"(.*\.?.*)\..*", s"\1") * "$suffix.$ext"


"""
    plot_histogram(bins, (f, f_style, f_name)...; kwargs...)

Plot histograms defined by their respective tuples `(f, f_style, f_name)` from histogram
values `f`, a PGFPlots style string `f_style`, and a legend entry `name`. The `bins` array
determines the positions on the x axis of each histogram.

**Keyword arguments**: `xlabel`, `ylabel`, `style`
"""
function plot_histogram{T1<:Number,T2<:Number}( bins   :: AbstractArray{T1,1},
                                                tups   :: Tuple{Array{T2,1},String,String}...;
                                                xlabel :: String = "levels",
                                                ylabel :: String = "edf",
                                                style  :: String = "",
                                                ymin   :: Float64 = 0.0 )
    # add ending value to complete last line segment
    levelwidth = bins[2] - bins[1]
    bins = vcat(bins[1]-1e-12, bins, bins[end] + levelwidth, bins[end] + levelwidth+1e-12) # add last level
    
    # map tuples to plots
    plots = map(tups) do tup
        f, f_style, f_name = tup # split tuple
        f = vcat(ymin, f, f[end], ymin) # ending value
        Plots.Linear(bins, f, legendentry = f_name, style = joinstyles("histogram plot", f_style))
    end
    
    # return axis object
    return Axis([plots...]; style = joinstyles("histogram axis", "ymode = log, log origin y=infty", style),
                xmin = bins[1], xmax = bins[end], xlabel = xlabel, ylabel = ylabel)
end

"""
    progress_style(i)

Obtain the style number `i` for progress plots.
"""
progress_style(i::Int) = @sprintf "progress plot, m%02d" i

"""
    joinstyles(styles...)

Join the given PGFPlots style strings.
"""
function joinstyles(styles::String...)
    a = IOBuffer()
    for s in styles
        length(s) > 0 && print(a, a.size > 0 ? ", " : "", s)
    end
    String(take!(a))
end


# res implementations
include("res/ibu.jl")
include("res/data.jl")
include("res/subsampling.jl")
include("res/smearing.jl")
include("res/clustering.jl")
include("res/weightfix.jl")
include("res/stepsize.jl")
include("res/smoothing.jl")
include("res/expand_reduce.jl")
include("res/comparison.jl")
include("res/time_series_contributions.jl")
include("res/smart_control.jl")


end

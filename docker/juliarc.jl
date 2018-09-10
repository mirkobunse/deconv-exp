# This file should contain site-specific commands to be executed on Julia startup
# Users should store their own personal commands in homedir(), in a file named .juliarc.jl

# initialize package repository
function install_dependencies()
    
    # install packages
    ENV["PYTHON"]=""
    map(Pkg.add, [ "Colors",
                   "DataFrames",
                   "Query",
                   "MLDataUtils",
                   "Discretizers",
                   "Distances",
                   "Distributions",
                   "HDF5",
                   "Images",
                   "Optim",
                   "PGFPlots",
                   "Polynomials",
                   "YAML",
                   "ScikitLearn",
                   "MultivariateStats" ])
    map(Pkg.clone, [ "git://github.com/mirkobunse/ComfyBase.jl.git",
                     "git://github.com/mirkobunse/ComfyCommons.jl.git",
                     "git://github.com/mirkobunse/CherenkovDeconvolution.jl.git" ])
    Pkg.checkout("CherenkovDeconvolution", "0.0.1", pull=false) # v0.0.1 is used in the experiments
    Pkg.build()
    
    # make directories
    for dir in ["data", "conf/job/gen", "res", "res/spectra"]
        try mkdir(joinpath(homedir(), "mt-exp", dir)) end
    end
    info("Volume successfully initialized. You should restart julia, now.")
    
end

# the usual
try
    using ComfyBase
    exit() = ComfyBase.confirmexit()
    ComfyBase.loadrc(confirm=isa(STDIN,Base.TTY)) # only ask when in a TTY
catch err
    warn("Volume is not yet initialized. You can do so by calling initvolume()")
end


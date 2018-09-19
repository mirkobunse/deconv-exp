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
                   "MultivariateStats" ])
    map(Pkg.clone, [ "https://github.com/mirkobunse/ComfyCommons.jl.git",
                     "https://github.com/mirkobunse/ScikitLearn.jl.git",
                     "https://github.com/mirkobunse/CherenkovDeconvolution.jl.git" ])
    Pkg.checkout("ComfyCommons", "julia-0.6", pull=false) # branch for julia's 0.6 version
    Pkg.checkout("ScikitLearn", "v0.4.0-fix", pull=false) # bug fix regarding pre-compilation
    Pkg.checkout("CherenkovDeconvolution", "0.0.1", pull=false) # v0.0.1 is used in the experiments
    Pkg.pin("DataFrames", v"0.11.6") # bug occurs in 0.11.7
    Pkg.build()
    
    # make directories
    for dir in ["data", "conf/job/gen", "res", "res/spectra"]
        try mkdir(joinpath(homedir(), "mt-exp", dir)) end
    end
    info("Volume successfully initialized. You should restart julia, now.")
    
end

# project-specific initialization
try
    isfile("_init.jl") && include(joinpath(pwd(), "_init.jl"))
catch err
    info("The volume is likely not initialized. Run `install_dependencies()` to solve this issue.")
end


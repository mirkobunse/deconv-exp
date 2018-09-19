# Pkg.clone("https://github.com/mirkobunse/ComfyCommons.jl.git")
# Pkg.checkout("ComfyCommons", "julia-0.6") # check out the deprecated branch for julia-0.6
using ComfyCommons

# Pkg.clone("https://github.com/mirkobunse/CherenkovDeconvolution.jl.git")
using ScikitLearn, CherenkovDeconvolution # docker container breaks if these are not imported, here
using ComfyCommons: info, warn
using DataFrames
using PGFPlots


MODULES = ["Util", "Data", "Conf", "Job", "Res"]
try ComfyCommons.Imports.importdir("src", modulenames=MODULES) end


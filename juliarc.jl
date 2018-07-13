using ComfyBase     # Pkg.clone("git://github.com/mirkobunse/ComfyBase.jl.git")
using ComfyCommons  # Pkg.clone("git://github.com/mirkobunse/ComfyCommons.jl.git")
using ScikitLearn, CherenkovDeconvolution # docker container breaks if these are not imported, here
using ComfyCommons: info, warn
using DataFrames
using PGFPlots


MODULES = ["Util", "Data", "Conf", "Job", "Res"]


try ComfyBase.importdir("src", modulenames=MODULES) end
Base.reload() = ComfyBase.reloadmodules(MODULES)


using JLD2
using Random
base_directory=pwd()
src_directory=joinpath(base_directory,"src")
scripts_dir=joinpath(base_directory,"scripts")
data_directory=joinpath(base_directory,"data")
readdir()

include(joinpath(src_directory,"FunctionalBayesExtremes.jl"))
using .FunctionalBayesExtremes


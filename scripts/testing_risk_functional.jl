using JLD2
using Random
base_directory=pwd()
src_directory=joinpath(base_directory,"src")
scripts_dir=joinpath(base_directory,"scripts")
data_directory=joinpath(base_directory,"data")


include(joinpath(src_directory,"FunctionalBayesExtremes.jl"))
using .FunctionalBayesExtremes


using Distributions

grid=default_Grid()
param=Parameter(α=2.0, β=0.5, c=3.0)
observation=Observation(param=param,grid=grid,num_runs=3000,num_sim=100)



N_est_c=100000
l2_a,l2_b=FunctionalBayesExtremes.l_2_fun_double_param(param_a=param, param_b=param, grid=grid, N_est_c=N_est_c, N_est_c_lim=N_est_c, quot_bound=10000000000.0, exceedance_observation=observation)
  


N_est_c=100001
l2_a,l2_b=FunctionalBayesExtremes.l_2_fun_double_param(param_a=param, param_b=param, grid=grid, N_est_c=N_est_c, N_est_c_lim=N_est_c, quot_bound=10000000000.0, exceedance_observation=observation)
       
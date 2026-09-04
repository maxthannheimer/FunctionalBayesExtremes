using JLD2
using Random
include("FunctionalBayesExtremes_MODULE.jl")
using .FunctionalBayesExtremes_MODULE
using Distributions

grid=default_Grid()
param=Parameter(α=2.0, β=0.5, c=3.0)
N_MCMC=10






for i in 1:1
#Generate observation
observation=Observation(param=param,grid=grid,num_runs=100,num_sim=100)
#Generate starting values via prior

start_alpha=gaussian_proposal(1.0,1.0)
start_beta=rand(Uniform(0.0, 2.0))
start_c=gaussian_proposal(1.0,1.5)






println("Time for approx: ")
@time (
dict_MCMC_approx=MCMC_approx_(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=Parameter(α=start_alpha, β=start_beta, c=start_c), grid=grid,N_est_c=400,N_cond_sim=60,N_est_d=60)
)

println("Time for single param: ")
@time (
dict_MCMC=MCMC_(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=Parameter(α=start_alpha, β=start_beta, c=start_c), grid=grid,N_est_c=800,N_cond_sim=120,N_est_d=120)
)

println("Time for double param: ")
@time (
dict_MCMC_double_param=MCMC_double_param(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=Parameter(α=start_alpha, β=start_beta, c=start_c), grid=grid,N_est_c=800,N_cond_sim=120,N_est_d=120)
)



save("TEST_Simulations_MCMC_single_param_"*randstring(20)*".jld2", "dict_MCMC",dict_MCMC,"dict_MCMC_approx",dict_MCMC_approx, "dict_MCMC_double_param",dict_MCMC_double_param, "grid",grid,"param",param,"observation",observation)
end



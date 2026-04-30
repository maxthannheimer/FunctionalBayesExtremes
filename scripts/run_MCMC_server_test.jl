using JLD2
using Random
include("FunctionalBayesExtremes_MODULE.jl")
using .FunctionalBayesExtremes_MODULE


grid=default_Grid()
param=Parameter(α=2.0, β=0.5, c=3.0)
N_MCMC=10

for i in 1:1
observation=Observation(param=param,grid=grid,num_runs=100,num_sim=100)
@time (
dict_MCMC=MCMC_(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=param, grid=grid,N_est_c=40000,N_cond_sim=100,N_est_d=600)
)

@time (
dict_MCMC_approx=MCMC_approx_(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=param, grid=grid,N_est_c=40000,N_cond_sim=100,N_est_d=600)
)

save("Simulations_MCMC_test_"*randstring(20)*".jld2", "dict_MCMC",dict_MCMC,"dict_MCMC_approx",dict_MCMC_approx,"grid",grid,"param",param,"observation",observation)
end



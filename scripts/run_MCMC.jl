using DrWatson
using JLD2
using Random
@quickactivate :FunctionalBayesExtremes

grid=small_default_Grid()
param=Parameter(α=0.5, β=1.5, c=3.0)
N_MCMC=10

observation=Observation(param=param,grid=grid,num_runs=1000,num_sim=100)

for i in 1:1
@time (
dict_MCMC=FunctionalBayesExtremes.MCMC_(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=param, grid=grid,N_est_c=20000,N_cond_sim=100,N_est_d=300)
)

@time (
dict_MCMC_approx=FunctionalBayesExtremes.MCMC_approx_(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=param, grid=grid,N_est_c=20000,N_cond_sim=100,N_est_d=300)
)

save("Simulations_MCMC_test"*randstring(20)*".jld2", "dict_MCMC",dict_MCMC,"dict_MCMC_approx",dict_MCMC_approx,"grid",grid,"param",param,"observation",observation)
end



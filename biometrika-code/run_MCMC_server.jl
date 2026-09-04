using JLD2
using Random
include("FunctionalBayesExtremes_MODULE.jl")
using .FunctionalBayesExtremes_MODULE
using Distributions

grid=default_Grid()
param=Parameter(α=2.0, β=0.5, c=3.0)
N_MCMC=10000



#Number of simulations running after each other
for i in 1:8
#Generate observation
observation=Observation(param=param,grid=grid,num_runs=5000,num_sim=100)
#Generate starting values via prior


start_alpha=gaussian_proposal(1.0,1.0)
start_beta=rand(Uniform(0.0, 2.0))
start_c=gaussian_proposal(1.0,1.5)


println("Params: ")
println("alpha: ", start_alpha," beta: ", start_beta, " c: ", start_c)

println("Time for approx: ")
@time (
dict_MCMC_approx=MCMC_approx_(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=Parameter(α=start_alpha, β=start_beta, c=start_c), grid=grid,N_est_c=40000,N_cond_sim=600,N_est_d=600)
)

#println("Time for single param: ")
#@time (
#dict_MCMC=MCMC_(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=Parameter(α=start_alpha, β=start_beta, c=start_c), grid=grid,N_est_c=80000,N_cond_sim=1200,N_est_d=1200)
#)

println("Time for double param: ")
@time (
dict_MCMC_double_param=MCMC_double_param(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=Parameter(α=start_alpha, β=start_beta, c=start_c), grid=grid,N_est_c=80000,N_cond_sim=1200,N_est_d=1200)
)



save("Simulations_MCMC_"*randstring(20)*".jld2", "dict_MCMC_approx",dict_MCMC_approx, "dict_MCMC_double_param",dict_MCMC_double_param, "grid",grid,"param",param,"observation",observation)
end



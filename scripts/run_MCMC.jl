using JLD2
using Random
include("/home/thannhmx/.julia/dev/FunctionalBayesExtremes/src/FunctionalBayesExtremes_MODULE.jl")
using .FunctionalBayesExtremes_MODULE
using Distributions

#Starting params

#Log Gaussian Code

c_start_prior=gaussian_proposal(1.0,1.5)
beta_start_prior=rand(Uniform(0.0, 2.0))
alpha_start_prior=gaussian_proposal(1.0,1.0)



grid=default_Grid()
param=Parameter(α=2.0, β=0.5, c=3.0)
N_MCMC=10
N_burn_in=1000
observation=Observation(param=param,grid=grid,num_runs=1000,num_sim=100)

#1000 Burn in
start_alpha=1.3829659530670813
start_beta=1.9899365928889978
start_c=1.317198054390836
#Parameter(α=start_alpha, β=start_beta, c=start_c)
param_start=Parameter(α=start_alpha, β=start_beta, c=start_c)
    num_obs=size(observation.obs_x0,1)
include("/home/thannhmx/.julia/dev/FunctionalBayesExtremes/src/functions.jl")
    param_vec = [Parameter(α=NaN, β=NaN, c=NaN) for i=1:N_MCMC+1]
    param_vec[1]=param_start
    number_exceed_vec = [NaN for i=1:N_MCMC]
    log_likelihood_vec = [NaN for i=1:N_MCMC]
    res_ell_X_vec = [[NaN for i in 1:num_obs] for j in 1:N_MCMC] 
    

trial=1
        (exceedance_observation, res_ell_X_vec[1]) = exceed_cond_sim_approx(observation=observation, threshold=1.0)

        #propose new params
        eps_beta=0.05 # 0.05 -> 0.025
        eps_c=0.1 #0.1 -> 0.05
        eps_alpha=0.1 #0.1 -> 0.05
 include("/home/thannhmx/.julia/dev/FunctionalBayesExtremes/src/MCMC.jl")       
        beta_eps,old_interval,new_interval=uniform_proposal(param_vec[trial].β,eps_beta,0.0,2.0)
        param_eps=Parameter(c=gaussian_proposal(param_vec[trial].c,eps_c),β=beta_eps,α=gaussian_proposal(param_vec[trial].α,eps_alpha))

      
        l1=l_1_fun(param=param_vec[trial], grid=grid,exceedance_observation=exceedance_observation)   
        l2=l_2_fun_approx(param=param_vec[trial],grid=grid,N_est_c=20000,exceedance_observation=exceedance_observation)
        l3=l_3_fun(exceedance_observation=exceedance_observation, param=param_vec[trial], threshold=1.0)
prior = log_likehood_log_gauss_1d_non_normalized(param_vec[trial].c,0.0,1.5)+log_likehood_log_gauss_1d_non_normalized(param_vec[trial].α,0.0,1.0)
        log_likelihood_old=sum([l1,l2,l3,prior,-log(old_interval)])





@time (
dict_MCMC_approx=MCMC_approx_(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=param_start, grid=grid,N_est_c=20000,N_cond_sim=100,N_est_d=300)


)


#for i in 1:1
@time (
dict_MCMC=MCMC_(N_MCMC=N_MCMC,observation=observation,threshold=1.0,param=param, grid=grid,N_est_c=20000,N_cond_sim=100,N_est_d=300)
)



save("TEST_Simulations"*randstring(20)*".jld2", "dict_MCMC",dict_MCMC,"dict_MCMC_approx",dict_MCMC_approx,"grid",grid,"param",param,"observation",observation)
#end

println("Time for Starting values: ")

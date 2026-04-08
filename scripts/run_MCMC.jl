using DrWatson
@quickactivate :FunctionalBayesExtremes

grid=default_Grid()
param=Parameter(α=1.5, β=1.5, c=4.0)
observation=Observation(param=param,grid=grid,num_runs=2000,num_sim=200)

dict_MCMC=FunctionalBayesExtremes.MCMC(N_MCMC=10,observation=observation,threshold=1.0,param=param, grid=grid,N_est_c=20000,N_cond_sim=1000,N_est_d=300)
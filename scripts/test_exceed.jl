using DrWatson
@quickactivate :FunctionalBayesExtremes

grid=default_Grid()
param=Parameter(α=1.5, β=1.5, c=4.0)
observation=Observation(param=param,grid=grid,num_runs=2000,num_sim=200)
size(observation.obs_data,1)
observation.obs_data[1,:]
test=FunctionalBayesExtremes.r_cond_gaussian(param=param,grid=grid,num_sim=10,observation=observation)

obs_exceed=FunctionalBayesExtremes.exceed_cond_sim(num_runs=10,observation=observation,threshold=1.0, param=param, grid=grid )
obs_exceed[1].obs_x0

typeof(obs_exceed)

FunctionalBayesExtremes.l_4_fun(observation=observation, threshold=1.0, param=param, grid=grid ,N_est_d=300)
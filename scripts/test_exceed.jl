using DrWatson
@quickactivate :FunctionalBayesExtremes

"""Here we check if the code for r_cond_W, exceed_cond_sim and the likelihood l_1_fun, l_2_fun, l_3_fun, l_4_fun functions runs without errors. """


"""Generate grid, parameters and observation data for testing if the code runs"""
grid=default_Grid()
param=Parameter(α=1.5, β=1.5, c=4.0)
observation=Observation(param=param,grid=grid,num_runs=2000,num_sim=200)

#test r_cond_W
test=FunctionalBayesExtremes.r_cond_gaussian_with_trafo(param=param,grid=grid,num_sim=10,observation=observation)
test2=FunctionalBayesExtremes.r_cond_W(param=param,grid=grid,num_sim=10,observation=observation)

obs_exceed=FunctionalBayesExtremes.exceed_cond_sim(num_runs=10,observation=observation,threshold=1.0, param=param, grid=grid )
obs_exceed[1].obs_x0

typeof(obs_exceed)

#test likelihood functions
FunctionalBayesExtremes.l_4_fun(exceedance_observation=observation, threshold=1.0, param=param, grid=grid ,N_est_d=300)
FunctionalBayesExtremes.l_4_fun(exceedance_observation=obs_exceed[1], threshold=1.0, param=param, grid=grid ,N_est_d=300)

FunctionalBayesExtremes.l_1_fun(param=param, grid=grid,exceedance_observation=observation)   
FunctionalBayesExtremes.l_1_fun(param=param, grid=grid,exceedance_observation=obs_exceed[1]) 

FunctionalBayesExtremes.l_3_fun(exceedance_observation=observation, param=param, threshold=1.0)
FunctionalBayesExtremes.l_3_fun(exceedance_observation=obs_exceed[1], param=param, threshold=1.0)

FunctionalBayesExtremes.l_2_fun(param=param,grid=grid,N_est_c=20000,exceedance_observation=obs_exceed[1])
FunctionalBayesExtremes.l_2_fun(param=param,grid=grid,N_est_c=20000,exceedance_observation=observation)


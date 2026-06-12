include("/home/thannhmx/.julia/dev/FunctionalBayesExtremes/src/FunctionalBayesExtremes.jl")
using .FunctionalBayesExtremes

#Import needed functions FBM_simu_fast and fBm for testing them

#define parameters and grid for testing
param_a=Parameter(α=1.0, β=1.9, c=3.0)
param_b=Parameter(α=1.0, β=1.5, c=3.0)
grid=default_Grid()
num_runs=100
num_sim=120
num_obs=100
#define observation for testing r_cond_gaussian_double_param
observation=Observation(param=param_a,grid=grid,num_runs=100,num_sim=num_obs)

gaussian_observation_a=FunctionalBayesExtremes.observation_trafo(observation=observation,param=param_a)
sigma_yy_inv_a = inv(FunctionalBayesExtremes.cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_coarse,  param=param_a, coord_x0=grid.coord_x0 )) #hier 
sigma_zy_a= FunctionalBayesExtremes.cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_fine, param=param_a, coord_x0=grid.coord_x0)'   
sigma_yy_inv_b = inv(FunctionalBayesExtremes.cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_coarse,  param=param_b, coord_x0=grid.coord_x0 )) #hier 
sigma_zy_b= FunctionalBayesExtremes.cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_fine, param=param_b, coord_x0=grid.coord_x0)'    


 res_a=[[[NaN for k in 1:grid.gridsize^2 ] for i in 1:num_sim ] for j in 1:size(gaussian_observation,1)]
 res_b=[[[NaN for k in 1:grid.gridsize^2 ] for i in 1:num_sim ] for j in 1:size(gaussian_observation,1)]

for j in 1:size(gaussian_observation,1)
    res_a[j],res_b[j]=FunctionalBayesExtremes.r_gaussian_double_param(param_a=param_a, param_b=param_b, grid=grid, num_sim=num_sim) 
end

res_old=[FunctionalBayesExtremes.r_gaussian(param=param_a, grid=grid, num_sim=num_sim) for j in 1:size(gaussian_observation_a,1)]

typeof(res_old) 
size(res_old,1) #num_sim
size(res_old[1],1) #num_obs
size(res_old[1][1],1) #gridsize

typeof(res_a) 
size(res_a,1) #num_sim
size(res_a[1],1) #num_obs
size(res_a[1][1],1) #gridsize

a,b=FunctionalBayesExtremes.r_cond_gaussian_double_param(param_a=param_a,param_b=param_b,grid=grid,num_sim=num_sim,gaussian_observation=gaussian_observation_a)
a
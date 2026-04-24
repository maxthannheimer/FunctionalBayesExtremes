""" Here different functions are defined, which are used in the rest of the code. """

""" variogram of the process 
 c⋅||x||^β """
function vario(;coord::Vector{Float64},param::Parameter)::Float64
    return param.c*sqrt(
        (coord[1])^2
        + (coord[2])^2
        )^param.β
end

""" vectorized version of the variogram function, for all grid points, normalized in x0. 
 c⋅||x-x0||^β for each location in coord_fine """

using LinearAlgebra

 function vec_vario_grid(;param::Parameter,grid::Grid)::Vector{Float64}
     param.c.*norm.(eachrow(grid.coord_fine.-grid.coord_x0')).^param.β
end

#TODO: test with tivial vectors by hand

""" vectorized version of the variogram function, for all grid points. 
 c⋅||x-x0||^β for each location in coord_vec, normalized in coord_x0 """
function vec_vario(;param::Parameter,coord_vec::Matrix{Float64},coord_x0::Vector{Float64})::Vector{Float64}
     param.c.*norm.(eachrow(coord_vec.-coord_x0')).^param.β
end


""" covariance function for two locations x and y 
 c⋅||x-x0||^β + c⋅||y-x0||^β -c⋅||x-y||^β , normalized in x0 """
function cov_fun_vario(;param::Parameter,coord_a::Vector{Float64},coord_b::Vector{Float64},coord_x0::Vector{Float64})::Float64
    return vario(coord=coord_a-coord_x0,param=param) + vario(coord=coord_b-coord_x0,param=param)-vario(coord=coord_a-coord_b,param=param)
end



"""calculate cov matrix for two matrices of coordinates and the variogramm given via param and normalized at ARBITRARY x0"""
function cov_mat_for_vectors(;coord_mat_a::Matrix{Float64}, coord_mat_b::Matrix{Float64}, param::Parameter, coord_x0::Vector{Float64})::Matrix{Float64}
    Nx=size(coord_mat_a,1)
    Ny=size(coord_mat_b,1)
    cov_mat=ones(Nx,Ny)
    for ix in 1:Nx
        for jy in 1:Ny
            cov_mat[ix,jy]=cov_fun_vario(param=param,coord_a=coord_mat_a[ix,:], coord_b=coord_mat_b[jy,:], coord_x0=coord_x0 )  
            #second implementation to check if everything works
            #cov_mat[ix,jy]=vario(coord_mat_a[ix,:]-coord_x0,param)+vario(coord_mat_b[jy,:]-coord_x0,param)-vario(coord_mat_a[ix,:]-coord_mat_b[jy,:],param)
        end
    end
    cov_mat
end


""" simulation of gaussian random vectors with brown resnick covariance, using the circulant embedding method.
#simulate numrep many exp(1/α[G(s)-G(x0)-γ(s-x0)])"""
function r_log_gaussian(;param::Parameter,grid::Grid,num_sim::Int) :: Vector{Vector{Float64}}
    cov_mat=cov_mat_for_vectors(coord_mat_a=grid.coord_fine, coord_mat_b=grid.coord_fine, param=param, coord_x0=grid.coord_x0)
    res = FBM_simu_fast(param=param, grid=grid, num_sim=num_sim)
    trend=vec_vario_grid(param=param,grid=grid)
    for i in 1:num_sim
            res[i] = exp.(1/param.α*(res[i] - trend .-res[i][grid.x0])) #variogram
    end
    res
end

function r_log_gaussian_alpha(;param::Parameter,grid::Grid,num_sim::Int) :: Vector{Vector{Float64}}
    cov_mat=cov_mat_for_vectors(coord_mat_a=grid.coord_fine, coord_mat_b=grid.coord_fine, param=param, coord_x0=grid.coord_x0)
    res = FBM_simu_fast(param=param, grid=grid, num_sim=num_sim)
    trend=vec_vario_grid(param=param,grid=grid)
    for i in 1:num_sim
            res[i] = exp.(1/param.α*(res[i] - trend .-res[i][grid.x0])) #variogram
    end
    res
end

""" simulation of gaussian random vectors with brown resnick covariance, using the circulant embedding method.
simulate numrep many [G(s)-G(x0)-γ(s-x0)] """
function r_gaussian(;param::Parameter,grid::Grid,num_sim::Int) :: Vector{Vector{Float64}}
    res = FBM_simu_fast(param=param, grid=grid, num_sim=num_sim)
    trend=vec_vario_grid(param=param,grid=grid)
    for i in 1:num_sim
            res[i] = res[i] - trend .-res[i][grid.x0]#variogram
    end
    res
end

""" simulation of log gaussian random vectors with brown resnick covariance, using the circulant embedding method.
simulate numrep many W(s)=exp(1/α [G(s)-G(x0)-γ(s-x0)] )"""
function r_W(;param::Parameter,grid::Grid,num_sim::Int) :: Vector{Vector{Float64}}
    res = FBM_simu_fast(param=param, grid=grid, num_sim=num_sim)
    trend=vec_vario_grid(param=param,grid=grid)
    for i in 1:num_sim
            res[i] = exp.(1/param.α*( res[i] - trend .-res[i][grid.x0])) #same as r_gaussian but with exp and 1/alpha
    end
    res
end


function r_gaussian_sparse(;param::Parameter,grid::Grid,num_sim::Int) :: Vector{Vector{Float64}}
    N=size(grid.coord_coarse,1)+1
    cov_matrix = zeros(Float64, N, N)
    coord_coarse_plus_x0 = vcat(grid.coord_coarse, grid.coord_x0')
    cov_matrix = cov_mat_for_vectors(coord_mat_a=coord_coarse_plus_x0, coord_mat_b=coord_coarse_plus_x0, param=param, coord_x0=grid.coord_x0).+1.0
    trend=vec_vario(param=param,coord_vec=coord_coarse_plus_x0,coord_x0=grid.coord_x0)
    
    # Generate a sample from the multivariate normal distribution
    mean_vector = zeros(Float64, N)
    #return cov_matrix
    res=[Vector{Float64}(undef,N) for i in 1:num_sim]
    for i in 1:num_sim
        res[i] = rand(MvNormal(mean_vector, cov_matrix))
    end
    for i in 1:num_sim
        res[i] = res[i] - trend .-res[i][end]   
    end
    return res
end


function r_W_sparse(;param::Parameter,grid::Grid,num_sim::Int) :: Vector{Vector{Float64}}
    N=size(grid.coord_coarse,1)+1
    cov_matrix = zeros(Float64, N, N)
    coord_coarse_plus_x0 = vcat(grid.coord_coarse, grid.coord_x0')
    cov_matrix = cov_mat_for_vectors(coord_mat_a=coord_coarse_plus_x0, coord_mat_b=coord_coarse_plus_x0, param=param, coord_x0=grid.coord_x0).+1.0
    trend=vec_vario(param=param,coord_vec=coord_coarse_plus_x0,coord_x0=grid.coord_x0)
    
    # Generate a sample from the multivariate normal distribution
    mean_vector = zeros(Float64, N)
    #return cov_matrix
    res=[Vector{Float64}(undef,N) for i in 1:num_sim]
    for i in 1:num_sim
        res[i] = rand(MvNormal(mean_vector, cov_matrix))
    end
    for i in 1:num_sim
        #res[i] = res[i] - trend .-res[i][end]  
        res[i] = exp.(1/param.α*( res[i] - trend .-res[i][end])) 
    end
    return res
end

""" conditional gaussian simulation (conditioning on arbitrary value vectors)"""

function r_cond_gaussian_observation_vectors(;param::Parameter,grid::Grid,num_sim::Int,cond_obs::Vector{Vector{Float64}}) #coord_x0 (hier c egal)
    #first dim number of simulated or observed data repetitions, second dim num_rep (how many simulations are wanted), third dim site in fine grid
    
    sigma_yy_inv = inv(cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_coarse,  param=param, coord_x0=grid.coord_x0 )) #hier 
    sigma_zy= cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_fine, param=param, coord_x0=grid.coord_x0)'   
    res=[r_gaussian(param=param, grid=grid, num_sim=num_sim) for j in 1:size(cond_obs,1)]
    for j in 1:size(cond_obs,1)
        for i in 1:num_sim
            normalized_coarse_observation = cond_obs[j][1:end-1].-cond_obs[j][end]
            res[j][i] =res[j][i] + sigma_zy*(sigma_yy_inv*(normalized_coarse_observation-res[j][i][grid.rows_coord_coarse])) #variogram
        end
    end
    res
end

""" conditional gaussian simulation (conditioning on transformed observations on coarse grid), result is still Gaussian"""
function r_cond_gaussian_with_trafo(;param::Parameter,grid::Grid,num_sim::Int,observation::Observation) :: Vector{Vector{Vector{Float64}}}
    #first dim number of simulated or observed data repetitions, second dim num_rep (how many simulations are wanted), third dim site in fine grid
    sigma_yy_inv = inv(cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_coarse,  param=param, coord_x0=grid.coord_x0 )) #hier 
    sigma_zy= cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_fine, param=param, coord_x0=grid.coord_x0)'   
    res=[r_gaussian(param=param, grid=grid, num_sim=num_sim) for j in 1:size(observation.obs_data,1)]
    for j in 1:size(observation.obs_data,1)
        for i in 1:num_sim
            #here observation is tranformed, sice we assume to observe W=exp(1/α G)=X/X_0 <=> G=α(log(X)-log(X_0)) 
            normalized_coarse_observation = param.α.*  (log.(observation.obs_data[j,:]).-log(observation.obs_x0[j]))
            res[j][i] =res[j][i] + sigma_zy*(sigma_yy_inv*(normalized_coarse_observation-res[j][i][grid.rows_coord_coarse])) #variogram
        end
    end
    res
end

""" conditional gaussian simulation (conditioning on transformed observations on coarse grid), result is still Gaussian"""
function r_cond_gaussian(;param::Parameter,grid::Grid,num_sim::Int,gaussian_observation::Vector{Vector{Float64}}) :: Vector{Vector{Vector{Float64}}}
    #first dim number of simulated or observed data repetitions, second dim num_rep (how many simulations are wanted), third dim site in fine grid
    sigma_yy_inv = inv(cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_coarse,  param=param, coord_x0=grid.coord_x0 )) #hier 
    sigma_zy= cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_fine, param=param, coord_x0=grid.coord_x0)'   
    res=[r_gaussian(param=param, grid=grid, num_sim=num_sim) for j in 1:size(gaussian_observation,1)]
        #grid.coord_fine,param,grid.coord_x0,num_sim,alpha) for j in 1:size(cond_obs,1)]
    for j in 1:size(gaussian_observation,1)
        for i in 1:num_sim
            res[j][i] =res[j][i] + sigma_zy*(sigma_yy_inv*(gaussian_observation[j]-res[j][i][grid.rows_coord_coarse])) #variogram
        end
    end
    res
end


function observation_trafo(;observation::Observation, param::Parameter)::Vector{Vector{Float64}}
    #transform observation to G space, since we observe W=exp(1/α G)=X/X_0 <=> G=α(log(X)-log(X_0)) 
    [param.α.*  (log.(observation.obs_data[j,:]).-log(observation.obs_x0[j])) for j in 1:size(observation.obs_data,1)]
end

#Simulation of W via conditional Gaussian simulation of G and transformation W=exp(1/α G)
function r_cond_W(;param::Parameter,grid::Grid,num_sim::Int,observation::Observation) :: Vector{Vector{Vector{Float64}}}
    tmp = r_cond_gaussian(param=param, grid=grid, num_sim=num_sim,gaussian_observation=observation_trafo(observation=observation,param=param))
    #Here Gaussian simulations of G are transformed to process W via W=exp(1/α G)
    return [[exp.(1/param.α.* w) for w in v]  for v in tmp] 
end

using Distributions

""" exceedance selection among all observations,
    output is Observation object containing only exceedances over threshold and risk functional value for each observation"""

function exceed_cond_sim(;num_runs::Int,observation::Observation,threshold::Float64, param::Parameter, grid::Grid )::Tuple{Observation, Vector{Float64}}
    num_obs=size(observation.obs_data,1)
    tmp_exp = r_cond_W(param=param, grid=grid, num_sim=num_runs+1,observation=observation)
    #Here Gaussian simulations of G are transformed to process W via W=exp(1/α G)
    res_ell_X = [0.0 for i in 1:num_obs] 

   #old_value=r_cond_log_gaussian(observation_data[1,:],observation_x0[1], coord_fine,coord_coarse,param,row_x0)
   for i in 1:num_obs #direkt num_obs viele simulations 
        old_value = tmp_exp[i][1]
        for trial in 1:num_runs
            proposal = tmp_exp[i][trial+1]
            acceptance_rate = min(1,mean(proposal)^param.α/mean(old_value)^param.α)   
            if (rand()< acceptance_rate)
                old_value=proposal
            end
        end
        res_ell_X[i]=observation.obs_x0[i]*mean(old_value)
    end
    if sum(res_ell_X.>threshold)==0
        Base._throw_argerror("not a single threshold exceedance")
    else
        #likelihood calculation and param updates
        #find all threshold exccedances and calculate the log of them
        ind=findall(res_ell_X.>threshold)
        exceed_obs=Observation(observation.sim_data, observation.obs_data[ind,:], observation.obs_x0[ind])
        (exceed_obs,res_ell_X)
    end
end

"""exceedance selection among all observations, but without the conditional simulation, only via the approximate risk functional, a.k.a. the mean of the coarse observations"""
function exceed_cond_sim_approx(;observation::Observation,threshold::Float64)::Tuple{Observation, Vector{Float64}}
    num_obs=size(observation.obs_data,1)
    res_ell_X = [0.0 for i in 1:num_obs] 
    for i in 1:num_obs #direkt num_obs viele simulations 
        res_ell_X[i]=mean(vcat(observation.obs_data[i,:],observation.obs_x0[i]))
    end
    if sum(res_ell_X.>threshold)==0
        Base._throw_argerror("not a single threshold exceedance")
    else
        #likelihood calculation and param updates
        #find all threshold exccedances and calculate the log of them
        ind=findall(res_ell_X.>threshold)
        exceed_obs=Observation(observation.sim_data, observation.obs_data[ind,:], observation.obs_x0[ind])
        (exceed_obs,res_ell_X)
    end
end

""" Now we build the likelihood parts according to our paper """

""" l4= ∏(for i in observation) P (x_i(s0)*r(W)>u | W=x_i/x_i(s0) ) 
    to estimate that emperically, we sample from W on the fine grid conditional on our normalized observations and then count how often r(W)>u/x(s0) holds"""

function l_4_fun(;exceedance_observation:: Observation, threshold::Float64, param::Parameter, grid::Grid ,N_est_d::Int)::Float64
    num_obs=size(exceedance_observation.obs_data,1)
    #Here simulations of W via W=exp(1/α G)
    tmp_exp = r_cond_W(param=param, grid=grid, num_sim=N_est_d,observation=exceedance_observation)
    res_est_prob = [0.0 for i in 1:num_obs] 
    for i in 1:num_obs #direkt num_obs viele simulations 
            counter = [mean(tmp_exp[i][trial]) for trial in 1:N_est_d].>threshold/exceedance_observation.obs_x0[i]
            res_est_prob[i] = sum(counter)/N_est_d
    end
    return sum(log.(res_est_prob))
end



""" l1= ∏(for i in observation) f_W(x_i/x_i(s0)) 
    use the exceedance observations and evaluate all the coarse gride densities f_W(x_i/x_i(s0)) for the log gaussian distrib W"""

""" this is the log gaussian density function for the coarse grid observations, which we use to calculate l1, the likelihood of the exceedance observations under the model W"""
function log_d_log_gaussian(mu,cov_mat_inv,x,inv_determinant)::Float64
   (-0.5*transpose(log.(x)-mu) * cov_mat_inv * (log.(x)-mu) )+0.5*log(inv_determinant)#+log(sqrt((2*pi)^(-size(x,1))))-sum(log.(x)) #hier sum log x, da wir die Dichte von W=exp(1/α G) berechnen, also müssen wir die Dichte von G mit der Ableitung von W nach G multiplizieren, was 1/W=1/exp(1/α G) ist, also log(1/exp(1/α G))=-log(W)=-log(x)
end

function l_1_fun(;param::Parameter, grid::Grid,exceedance_observation::Observation)::Float64
    cov_mat_coarse_inv= param.α^2*inv(cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_coarse, param=param, coord_x0=grid.coord_x0)) #hier 
    inv_determinant = det(cov_mat_coarse_inv)
    trend = -vec_vario(param=param,coord_vec=grid.coord_coarse,coord_x0=grid.coord_x0)/param.α
    sum([(log_d_log_gaussian(trend, cov_mat_coarse_inv, exceedance_observation.obs_data[i,:]./exceedance_observation.obs_x0[i], inv_determinant)) for i in 1:size(exceedance_observation.obs_data,1)])
end



""" l3 = ∏( for i in observation) α (x_i(s_0)/thresh)^(-α-1) """

function l_3_fun(;exceedance_observation::Observation, param::Parameter, threshold::Float64)::Float64
    sum([log(param.α)+log((exceedance_observation.obs_x0[i]/threshold))*(-param.α-1) for i in 1:size(exceedance_observation.obs_x0,1)])
end

""" l2= est(1/C) * size(obs) """

function l_2_fun(;param::Parameter, grid::Grid, N_est_c::Int, exceedance_observation::Observation)::Float64
    tmp = r_W(param = param, grid = grid, num_sim = N_est_c) 
    -size(exceedance_observation.obs_x0,1) * log(mean([mean(tmp[i] 
        )^(param.α) for i in 1:N_est_c]))  
     # minus for 1/c_l (in log)
 end

 """ l2_approx= est_approx(1/C) * size(obs) """
function l_2_fun_approx(;param::Parameter, grid::Grid, N_est_c::Int, exceedance_observation::Observation)::Float64
    tmp = r_W_sparse(param = param, grid = grid, num_sim = N_est_c) 
    -size(exceedance_observation.obs_x0,1) * log(mean([mean(tmp[i] 
        )^(param.α) for i in 1:N_est_c]))  
     # minus for 1/c_l (in log)
 end

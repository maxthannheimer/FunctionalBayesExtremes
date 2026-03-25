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

#param=Parameter(α=1.0, β=2.0, c=3.0)
#grid=Grid()
#grid.gridsize
#vec_vario_grid(param=param,grid=grid)
#vec_vario(param=param,coord_vec=grid.coord_fine,coord_x0=grid.coord_x0)

""" covariance function for two locations x and y 
 c⋅||x-x0||^β + c⋅||y-x0||^β -c⋅||x-y||^β , normalized in x0 """
function cov_fun_vario(;param::Parameter,coord_a::Vector{Float64},coord_b::Vector{Float64},coord_x0::Vector{Float64})::Float64
    return vario(coord=coord_a-coord_x0,param=param) + vario(coord=coord_b-coord_x0,param=param)-vario(coord=coord_a-coord_b,param=param)
end

#cov_fun_vario(param=param,coord_a=grid.coord_fine[5,:],coord_b=grid.coord_fine[2,:],coord_x0=grid.coord_x0)


"""calculate cov matrix for two matrices of coordinates and the variogramm given via param and normalized at x0"""
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

#cov_fun_vario(param=param,coord_a=grid.coord_fine[1,:], coord_b=grid.coord_fine[2,:], coord_x0=grid.coord_x0 )  
#cov_mat_for_vectors(coord_mat_a=grid.coord_fine, coord_mat_b=grid.coord_fine, param=param, coord_x0=grid.coord_x0)

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

""" simulation of gaussian random vectors with brown resnick covariance, using the circulant embedding method.
simulate numrep many 1/α[G(s)-G(x0)-γ(s-x0)] """
function r_gaussian(;param::Parameter,grid::Grid,num_sim::Int) :: Vector{Vector{Float64}}
    res = FBM_simu_fast(param=param, grid=grid, num_sim=num_sim)
    trend=vec_vario_grid(param=param,grid=grid)
    for i in 1:num_sim
            res[i] = 1/param.α*(res[i] - trend .-res[i][grid.x0])#variogram
    end
    res
end

function r_gaussian_sparse(;param::Parameter,coord_coarse::Matrix{Float64},coord_x0::Vector{Float64},num_sim::Int) :: Vector{Vector{Float64}}
    N=size(coord_coarse,1)+1
    cov_matrix = zeros(Float64, N, N)
    coord_coarse_plus_x0 = vcat(coord_coarse, coord_x0')
    #for i in 1:N
     #   for j in 1:N
      #      cov_matrix[i, j] =  param.c*((norm(coord_coarse_plus_x0[i,:])^param.β + norm(coord_coarse_plus_x0[j,:])^param.β - norm(coord_coarse_plus_x0[i,:] - coord_coarse_plus_x0[j,:])^param.β)).+1.0
    #    end
    #end
    cov_matrix = cov_mat_for_vectors(coord_mat_a=coord_coarse_plus_x0, coord_mat_b=coord_coarse_plus_x0, param=param, coord_x0=coord_x0).+1.0
    trend=vec_vario(param=param,coord_vec=coord_coarse_plus_x0,coord_x0=coord_x0)
    
    # Generate a sample from the multivariate normal distribution
    mean_vector = zeros(Float64, N)
    #return cov_matrix
    res=[Vector{Float64}(undef,N) for i in 1:num_sim]
    for i in 1:num_sim
        res[i] = rand(MvNormal(mean_vector, cov_matrix))
    end
    for i in 1:num_sim
            res[i] = 1/param.α*(res[i] - trend .-res[i][end])#variogram
    end
    #  res2=[Vector{Float64}(undef,N) for i in 1:num_sim]
    # for i in 1:num_sim
    #     res2[i] = rand(MvNormal(mean_vector, cov_matrix))
    # end
    # for i in 1:num_sim
    #         res2[i] = 1/param.α*(res2[i] - trend .-res2[i][end])#variogram
    # end
    return res#cov_matrix,cov_matrix_2
end



""" conditional gaussian simulation (conditioning on the coarse grid observations)"""
function r_cond_gaussian(;param::Parameter,grid::Grid,num_sim::Int,cond_obs::Vector{Vector{Float64}}) #coord_x0 (hier c egal)
    #first dim number of simulated or observed data repetitions, second dim num_rep (how many simulations are wanted), third dim site in fine grid
    
    sigma_yy_inv = inv(cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_coarse,  param=param, coord_x0=grid.coord_x0 )) #hier 
    sigma_zy= cov_mat_for_vectors(coord_mat_a=grid.coord_coarse, coord_mat_b=grid.coord_fine, param=param, coord_x0=grid.coord_x0)'   
    res=[r_gaussian(param=param, grid=grid, num_sim=num_sim) for j in 1:size(cond_obs,1)]
        #grid.coord_fine,param,grid.coord_x0,num_sim,alpha) for j in 1:size(cond_obs,1)]
    for j in 1:size(cond_obs,1)
        for i in 1:num_sim
            normalized_coarse_observation = cond_obs[j][1:end-1].-cond_obs[j][end]
            res[j][i] =res[j][i] + sigma_zy*(sigma_yy_inv*(normalized_coarse_observation-res[j][i][grid.rows_coord_coarse])) #variogram
        end
    end
    res
end

#TODO: simulate conditional on random data and see if the empirical covariance is identical
#add observation and use it for cond sim
#grid.coord_fine
#grid.coord_coarse
#grid.coord_x0
#grid.x0

#grid_test=Grid(gridsize=10, N_coarse=5)
#grid_test.rows_coord_coarse
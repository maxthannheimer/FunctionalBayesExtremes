using DrWatson
using Distributions
using LinearAlgebra
@quickactivate :FunctionalBayesExtremes

""" Generate some observation data and save it to the data directory. """



grid=default_Grid()
param=Parameter(α=1.5, β=1.5, c=4.0)
# i,d= first(enumerate(dicts))
# param=Parameter(α=d["α"], β=d["β"], c=d["c"])
#cov_matrix_2 = cov_mat_for_vectors(coord_mat_a=coord_coarse_plus_x0, coord_mat_b=coord_coarse_plus_x0, param=param, coord_x0=grid.coord_x0).+1.0
#cholesky(cov_matrix_2)
#eigen(cov_matrix_2)


num_sim=20000
@time (coarse_gauss_data=FunctionalBayesExtremes.r_gaussian_sparse(param=param,coord_coarse=grid.coord_coarse,coord_x0=grid.coord_x0,num_sim=num_sim))




@time (fine_gauss_data=FunctionalBayesExtremes.r_cond_gaussian_observation_vectors(param=param,grid=grid,num_sim=1,cond_obs=coarse_gauss_data))


#standard method using cov matrix
#FBM_covmat_res=fBm(param=param,grid=grid,num_sim=num_sim)
emp_cov_mat=zeros(grid.gridsize^2,grid.gridsize^2)
mean_vec=zeros(grid.gridsize^2)

#calculate mean STANDARD METHOD
for i in 1:grid.gridsize^2
    mean_vec[i]=mean([fine_gauss_data[rep][1][i] for rep in 1:num_sim])
end

#calculate empirical covariance matrix STANDARD METHOD
for row in 1:grid.gridsize^2
    for col in 1:grid.gridsize^2
        emp_cov_mat[row,col]=1/(num_sim-1)*sum([(fine_gauss_data[rep][1][row]-mean_vec[row])*(fine_gauss_data[rep][1][col]-mean_vec[col]) for rep in 1:num_sim])
    end
end



#calculate true covariance matrix
true_cov_mat=zeros(grid.gridsize^2,grid.gridsize^2)
for row in 1:grid.gridsize^2
    for col in 1:grid.gridsize^2
        true_cov_mat[row,col]=param.c*(norm(grid.coord_fine[row,:]-grid.coord_x0)^param.β+norm(grid.coord_fine[col,:]-grid.coord_x0)^param.β-norm(grid.coord_fine[row,:]-grid.coord_fine[col,:])^param.β)
    end
end

#safe and load results
safesave(datadir("cov_mat_simulations_cond_sim", savename("cond_sim_empircical_covmat",num_sim,"jld2")), "emp_cov_mat",emp_cov_mat,"true_cov_mat",true_cov_mat,"param",param,"grid",grid,"fine_gauss_data",fine_gauss_data,"coarse_gauss_data",coarse_gauss_data)
emp_cov_mat, true_cov_mat, param, grid,fine_gauss_data,coarse_gauss_data=load(datadir("cov_mat_simulations_cond_sim", savename("cond_sim_empircical_covmat",num_sim,"jld2")), "emp_cov_mat", "true_cov_mat", "param", "grid","fine_gauss_data","coarse_gauss_data")





#||.||_∞ error between empirical and true covariance matrix: Conditional Gaussian Simulation
maximum(
    abs.((emp_cov_mat-true_cov_mat)  )#./(true_cov_mat.+10^(-8)))
    )


true_cov_mat2=FunctionalBayesExtremes.cov_mat_for_vectors(;coord_mat_a=grid.coord_fine, coord_mat_b=grid.coord_fine, param=param, coord_x0=grid.coord_x0)

maximum(
    abs.(emp_cov_mat-true_cov_mat2)
)



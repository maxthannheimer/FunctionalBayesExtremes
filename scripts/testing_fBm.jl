""" Testing fractional Brownian motion simulation by calculating the empirical covariance matrix and comparing it to the true covariance matrix. """

using DrWatson
@quickactivate :FunctionalBayesExtremes

using .FunctionalBayesExtremes

#Import needed functions FBM_simu_fast and fBm for testing them

#define parameters and grid for testing
param=Parameter(α=1.0, β=1.9, c=3.0)
grid=Grid()


#simulate num_sim many fBm samples and calculate empirical covariance matrix, compare to true covariance matrix
num_sim=100

#circulant embedding method
FBM_res=FunctionalBayesExtremes.FBM_simu_fast(param=param,grid=grid,num_sim=num_sim)
emp_cov_mat=zeros(grid.gridsize^2,grid.gridsize^2)
mean_vec=zeros(grid.gridsize^2)

#calculate mean CIRCULANT EMBEDDING
for i in 1:grid.gridsize^2
    mean_vec[i]=mean([FBM_res[rep][i] for rep in 1:num_sim])
end
#calculate empirical covariance matrix CIRCULANT EMBEDDING
for row in 1:grid.gridsize^2
    for col in 1:grid.gridsize^2
        emp_cov_mat[row,col]=1/(num_sim-1)*sum([(FBM_res[rep][row]-mean_vec[row])*(FBM_res[rep][col]-mean_vec[col]) for rep in 1:num_sim])
    end
end



#standard method using cov matrix
FBM_covmat_res=fBm(param=param,grid=grid,num_sim=num_sim)
emp_cov_mat2=zeros(grid.gridsize^2,grid.gridsize^2)
mean_vec2=zeros(grid.gridsize^2)

#calculate mean STANDARD METHOD
for i in 1:grid.gridsize^2
    mean_vec2[i]=mean([FBM_covmat_res[rep][i] for rep in 1:num_sim])
end

#calculate empirical covariance matrix STANDARD METHOD
for row in 1:grid.gridsize^2
    for col in 1:grid.gridsize^2
        emp_cov_mat2[row,col]=1/(num_sim-1)*sum([(FBM_covmat_res[rep][row]-mean_vec2[row])*(FBM_covmat_res[rep][col]-mean_vec2[col]) for rep in 1:num_sim])
    end
end



#calculate true covariance matrix
true_cov_mat=zeros(grid.gridsize^2,grid.gridsize^2)
for row in 1:grid.gridsize^2
    for col in 1:grid.gridsize^2
        true_cov_mat[row,col]=param.c*(norm(grid.coord_fine[row,:])^param.β+norm(grid.coord_fine[col,:])^param.β-norm(grid.coord_fine[row,:]-grid.coord_fine[col,:])^param.β)
    end
end






#safe and load results
safesave(datadir("cov_mat_simulations", savename("FBM_empircical_covmat",num_sim,"jld2")), "emp_cov_mat",emp_cov_mat,"emp_cov_mat2",emp_cov_mat2,"true_cov_mat",true_cov_mat,"param",param,"grid",grid)
load(datadir("cov_mat_simulations", savename("FBM_empircical_covmat",num_sim,"jld2")), "emp_cov_mat", "emp_cov_mat2", "true_cov_mat", "param", "grid")





#||.||_∞ error between empirical and true covariance matrix: CIRCULANT EMBEDDING
maximum(
    abs.((emp_cov_mat-true_cov_mat)  )#./(true_cov_mat.+10^(-8)))
    )

#||.||_∞ error between empirical and true covariance matrix: CIRCULANT EMBEDDING
maximum(
    abs.((emp_cov_mat2-true_cov_mat) )#./(true_cov_mat.+10^(-8)))
    )

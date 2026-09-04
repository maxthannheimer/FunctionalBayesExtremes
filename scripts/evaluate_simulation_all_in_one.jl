using DrWatson
using Distributions
#@quickactivate :FunctionalBayesExtremes
include(srcdir("FunctionalBayesExtremes_MODULE.jl"))
using .FunctionalBayesExtremes_MODULE
using PrettyTables
""" RSME, BIAS, empirical coverage and interval width calculation for the conditional simulation and the approximation method"""

####################
#set name of the simulation results folder
#date_string="2026_05_06"
#date_string="2026_06_08"
#date_string="2026_08_07"
date_string="2026_07_simulation_study"
#date_string="2026_06_15_single_param"
####################
N_burn_in=1000
quantile_val=0.1
#true_param=Parameter(α=0.5, β=1.5, c=3.0)
true_param=Parameter(α=2.0, β=0.5, c=3.0)
true_param_dict=Dict("α" => true_param.α, "β" => true_param.β, "c" => true_param.c)
total_simulation_number=size(readdir(datadir("exp_raw",date_string)),1)



#empty dicts to store results and evaluation metrics
est_cond_sim_mean=Dict("β" => [NaN for i in 1:total_simulation_number], 
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_cond_sim_sd=Dict("β" => [NaN for i in 1:total_simulation_number], 
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_cond_sim_median=Dict( "β" => [NaN for i in 1:total_simulation_number], 
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_cond_sim_lower_quantile=Dict( "β" => [NaN for i in 1:total_simulation_number], 
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_cond_sim_upper_quantile=Dict( "β" => [NaN for i in 1:total_simulation_number], 
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_cond_sim_interval_width=Dict( "β" => [NaN for i in 1:total_simulation_number], 
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_cond_sim_empirical_coverage=Dict( "β" => [NaN for i in 1:total_simulation_number], 
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_cond_sim_normal_coverage=Dict( "β" => [NaN for i in 1:total_simulation_number], 
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])


    RMSE_cond_sim_mean=Dict( "β" => NaN, 
                        "c" => NaN, 
                        "α" => NaN)
    RMSE_cond_sim_median=Dict( "β" => NaN,
                        "c" => NaN, 
                        "α" => NaN)
    BIAS_cond_sim_mean=Dict( "β" => NaN, 
                        "c" => NaN, 
                        "α" => NaN)
    BIAS_cond_sim_median=Dict( "β" => NaN, 
                        "c" => NaN, 
                        "α" => NaN)
    RMSE_cond_sim_median=Dict( "β" => NaN,
                        "c" => NaN, 
                        "α" => NaN)
    cond_sim_interval_width=Dict( "β" => NaN,
                        "c" => NaN, 
                        "α" => NaN)
    cond_sim_empirical_coverage=Dict( "β" => NaN,
                        "c" => NaN, 
                        "α" => NaN)
    cond_sim_normal_coverage=Dict( "β" => NaN,
                        "c" => NaN, 
                        "α" => NaN)


pretty_table_data=Dict{}()
#sim_data_names=["dict_MCMC_approx", "dict_MCMC","dict_MCMC_double_param"]
sim_data_names=["dict_MCMC_approx","dict_MCMC_double_param"]

for sim_data_name in sim_data_names                    

#here starts the actual calculation
for file_number in 1:size(readdir(datadir("exp_raw",date_string)),1)
#file_number=1 
file_string=readdir(datadir("exp_raw",date_string))
    data_dict_tmp=load(datadir("exp_raw", date_string, file_string[file_number]))
    param_res_dict=Dict(
        "α" => [param.α for param in data_dict_tmp[sim_data_name]["param"]],
        "β" => [param.β for param in data_dict_tmp[sim_data_name]["param"]],
        "c" => [param.c for param in data_dict_tmp[sim_data_name]["param"]]
    )
    for key in keys(param_res_dict)          
                    est_cond_sim_mean[key][file_number]=mean(param_res_dict[key][N_burn_in:end])
                    est_cond_sim_sd[key][file_number]=std(param_res_dict[key][N_burn_in:end])
                    est_cond_sim_median[key][file_number]=median(param_res_dict[key][N_burn_in:end])
                    est_cond_sim_lower_quantile[key][file_number]=quantile(param_res_dict[key][N_burn_in:end], quantile_val)
                    est_cond_sim_upper_quantile[key][file_number]=quantile(param_res_dict[key][N_burn_in:end], 1-quantile_val)

                    est_cond_sim_empirical_coverage[key][file_number]= ((est_cond_sim_lower_quantile[key][file_number] .<= true_param_dict[key] ) .* (est_cond_sim_upper_quantile[key][file_number] .>= true_param_dict[key] ))*1.0
                    #est_cond_sim_empirical_coverage[key][number_parallel_int,number_sim_int]=sum( (est_cond_sim_lower_quantile[key][number_parallel_int,number_sim_int] .<= true_param[key] ) .* (est_cond_sim_upper_quantile[key][number_parallel_int,number_sim_int] .>= true_param[key] ) )*1.0
                  
                    (lower_quantile, upper_quantile)=quantile.(Normal(est_cond_sim_mean[key][file_number],est_cond_sim_sd[key][file_number]), [quantile_val, 1-quantile_val])
                    est_cond_sim_normal_coverage[key][file_number]=( lower_quantile .<= true_param_dict[key] ) .* (upper_quantile .>= true_param_dict[key] )*1.0
                   est_cond_sim_interval_width[key][file_number]=est_cond_sim_upper_quantile[key][file_number]-est_cond_sim_lower_quantile[key][file_number]
 end

    for key in keys(est_cond_sim_mean)
                RMSE_cond_sim_mean[key]=round(mean((est_cond_sim_mean[key].-true_param_dict[key]).^2)^0.5 , digits=3)
                RMSE_cond_sim_median[key]=round(mean((est_cond_sim_median[key].-true_param_dict[key]).^2)^0.5 , digits=3)
                BIAS_cond_sim_mean[key]=round(mean(est_cond_sim_mean[key].-true_param_dict[key]) , digits=3)
                BIAS_cond_sim_median[key]=round(mean(est_cond_sim_median[key].-true_param_dict[key]), digits=3)

                cond_sim_empirical_coverage[key]=round(mean(est_cond_sim_empirical_coverage[key]), digits=3)
                cond_sim_normal_coverage[key]=round(mean(est_cond_sim_normal_coverage[key]), digits=3)
                cond_sim_interval_width[key]=round(mean((est_cond_sim_interval_width[key]).^2)^0.5 , digits=3)
            
    end
end




   
""" print results in a pretty table"""   
key_array=["β", "c", "α"]
    for key in key_array
        pretty_table_data[key*sim_data_name]=[sim_data_name,key,true_param_dict[key],RMSE_cond_sim_mean[key], RMSE_cond_sim_median[key], BIAS_cond_sim_mean[key], BIAS_cond_sim_median[key], cond_sim_empirical_coverage[key], cond_sim_normal_coverage[key]]
    end

end
keys(pretty_table_data)
data=permutedims(hcat([pretty_table_data[key] for key in keys(pretty_table_data)]...))
data=data[sortperm(data[:, 3]), :] 
column_labels =[["simulation type","parameter","Parameter_value", "RMSE_Mean", "RMSE_Median", "BIAS_Mean", "BIAS_Median", "Emp_Coverage", "Normal_Coverage"]]
println(date_string)
pretty_table(data; column_labels)





typeof(key_array[1])














date_string="2026_06_15_double_param"
date_string="2026_06_15_single_param"


mean_vec=[NaN for i in 1:size(readdir(datadir("exp_raw",date_string)),1)]
mean_vec_approx=[NaN for i in 1:size(readdir(datadir("exp_raw",date_string)),1)]

max_vec=[NaN for i in 1:size(readdir(datadir("exp_raw",date_string)),1)]
max_vec_approx=[NaN for i in 1:size(readdir(datadir("exp_raw",date_string)),1)]
min_vec=[NaN for i in 1:size(readdir(datadir("exp_raw",date_string)),1)]
min_vec_approx=[NaN for i in 1:size(readdir(datadir("exp_raw",date_string)),1)]
#here starts the actual calculation
for file_number in 1:size(readdir(datadir("exp_raw",date_string)),1)
#file_number=1 
file_string=readdir(datadir("exp_raw",date_string))
    data_dict_tmp=load(datadir("exp_raw", date_string, file_string[file_number]))
    param_res_dict=Dict(
        "α" => [param.α for param in data_dict_tmp["dict_MCMC"]["param"]],
        "β" => [param.β for param in data_dict_tmp["dict_MCMC"]["param"]],
        "c" => [param.c for param in data_dict_tmp["dict_MCMC"]["param"]]
    )
    param_res_dict_approx=Dict(
        "α" => [param.α for param in data_dict_tmp["dict_MCMC_approx"]["param"]],
        "β" => [param.β for param in data_dict_tmp["dict_MCMC_approx"]["param"]],
        "c" => [param.c for param in data_dict_tmp["dict_MCMC_approx"]["param"]]
    )
param_string="c"


#println("file_number: ", file_number, ", mean: ", mean(param_res_dict[param_string]),", min: ", minimum(param_res_dict[param_string]),", mean: ",mean(param_res_dict[param_string]))
#println("file_number: ", file_number, ", mean: ", mean(param_res_dict_approx[param_string]),", min: ", minimum(param_res_dict_approx[param_string]),", mean: ",mean(param_res_dict_approx[param_string]))
max_vec[file_number]=maximum(param_res_dict[param_string])
max_vec_approx[file_number]=maximum(param_res_dict_approx[param_string])
mean_vec[file_number]=mean(param_res_dict[param_string])
mean_vec_approx[file_number]=mean(param_res_dict_approx[param_string])
min_vec[file_number]=minimum(param_res_dict[param_string])
min_vec_approx[file_number]=minimum(param_res_dict_approx[param_string])

end

using Plots
histogram(max_vec, title="Max of c "*date_string, label="MCMC", xlabel="Maximum of c", ylabel="Frequency")
histogram(max_vec_approx, title="Max c (approx)"*date_string, label="MCMC_approx", xlabel="Maximum of c", ylabel="Frequency")
histogram(mean_vec, title="Mean of c"*date_string, label="MCMC", xlabel="Mean of c", ylabel="Frequency")
histogram(mean_vec_approx, title="Mean of c (approx)"*date_string, label="MCMC_approx", xlabel="Mean of c", ylabel="Frequency")


histogram(min_vec, title="Min of c "*date_string, label="MCMC", xlabel="Minimum of c", ylabel="Frequency")
histogram(min_vec_approx, title="Min c (approx)"*date_string, label="MCMC_approx", xlabel="Minimum of c", ylabel="Frequency")

#mean_vec_double_param=mean_vec
#mean_vec_approx_double_param=mean_vec_approx


#mean_vec_single_param=mean_vec
#mean_vec_approx_single_param=mean_vec_approx

sort(mean_vec_double_param)
sort(mean_vec_single_param)

"""Plotting α,β and c for different realisations of the MCMC and MCMC_approx chains"""


file_number=0



    file_number=file_number+1
    file_number=29
    file_string=readdir(datadir("exp_raw", date_string))
    data_dict_tmp=load(datadir("exp_raw", date_string, file_string[file_number]))
    param_res_dict=Dict(
        "α" => [param.α for param in data_dict_tmp["dict_MCMC"]["param"]],
        "β" => [param.β for param in data_dict_tmp["dict_MCMC"]["param"]],
        "c" => [param.c for param in data_dict_tmp["dict_MCMC"]["param"]]
    )
    param_res_dict_approx=Dict(
        "α" => [param.α for param in data_dict_tmp["dict_MCMC_approx"]["param"]],
        "β" => [param.β for param in data_dict_tmp["dict_MCMC_approx"]["param"]],
        "c" => [param.c for param in data_dict_tmp["dict_MCMC_approx"]["param"]]
    )

    plots = Vector{}(undef, 6)
    plots[1] = scatter(1:length(param_res_dict["c"]), param_res_dict["c"], title="MCMC", label="c");
    hline!([mean(param_res_dict["c"][N_burn_in:end])], label="MCMC mean", color=:red);
    plots[2] = scatter(1:length(param_res_dict["β"]), param_res_dict["β"], title="MCMC", label="β");
    hline!([mean(param_res_dict["β"][N_burn_in:end])], label="MCMC mean", color=:red);
    plots[3] = scatter(1:length(param_res_dict["α"]), param_res_dict["α"], title="MCMC", label="α");
    hline!([mean(param_res_dict["α"][N_burn_in:end])], label="MCMC mean", color=:red);
    plots[4] = scatter(1:length(param_res_dict_approx["c"]), param_res_dict_approx["c"], title="MCMC_approx", label="c");
    hline!([mean(param_res_dict_approx["c"][N_burn_in:end])], label="MCMC_approx mean", color=:red);
    plots[5] = scatter(1:length(param_res_dict_approx["β"]), param_res_dict_approx["β"], title="MCMC_approx", label="β");
    hline!([mean(param_res_dict_approx["β"][N_burn_in:end])], label="MCMC_approx mean", color=:red);
    plots[6] = scatter(1:length(param_res_dict_approx["α"]), param_res_dict_approx["α"], title="MCMC_approx", label="α");
    hline!([mean(param_res_dict_approx["α"][N_burn_in:end])], label="MCMC_approx mean", color=:red);
    

    plot(plots[1],plots[4],plots[2],plots[5],plots[3],plots[6], layout=(3,2),size=(1200, 800))

    #Plot Number of Exceed and likelhood
    plot(1:length(data_dict_tmp["dict_MCMC"]["Number of exceedance"]), data_dict_tmp["dict_MCMC"]["Number of exceedance"])
    plot(1:length(data_dict_tmp["dict_MCMC"]["log_likelihood"]), data_dict_tmp["dict_MCMC"]["log_likelihood"])




#check for large deviations in c estimation
sort(est_cond_sim_mean["c"])
sort(est_approx_mean["c"])

for i in 1:size(readdir(datadir("exp_raw", date_string)),1)
if est_cond_sim_mean["c"][i] > 4
    println("Simulation number: ", i, " has a mean estimate of c larger than 5: ", est_cond_sim_mean["c"][i])
end
end



# β   
    println("β")
    mean(param_res_dict["β"][N_burn_in:end])
    mean(param_res_dict_approx["β"][N_burn_in:end])
# α
    println("α")
    mean(param_res_dict["α"][N_burn_in:end])
    mean(param_res_dict_approx["α"][N_burn_in:end])
# c
    println("c")
    mean(param_res_dict["c"][N_burn_in:end])
    mean(param_res_dict_approx["c"][N_burn_in:end])




"""Repeated l2 estimation for approx and cond sim, then histogram, mean and standarddeviation"""


observation=data_dict_tmp["observation"]
grid=data_dict_tmp["grid"]  
N_hist=300
coarse_est_vec=   [    FunctionalBayesExtremes.l_2_fun_approx(param=true_param,grid=grid,N_est_c=20000,exceedance_observation=observation) for i in 1:N_hist]

@time (
fine_est_vec=   [    FunctionalBayesExtremes.l_2_fun(param=true_param,grid=grid,N_est_c=20000,exceedance_observation=observation) for i in 1:N_hist]
)


save(datadir("exp_processed", date_string, "coarse_est_vec.jld2"), "coarse_est_vec", coarse_est_vec)
save(datadir("exp_processed", date_string, "fine_est_vec.jld2"), "fine_est_vec", fine_est_vec)
coarse_est_vec=load(datadir("exp_processed", date_string, "coarse_est_vec.jld2"), "coarse_est_vec")
fine_est_vec=load(datadir("exp_processed", date_string, "fine_est_vec.jld2"), "fine_est_vec")


histogram(coarse_est_vec, title="Coarse estimation of l_2", label="l_2 coarse estimates", xlabel="l_2 estimate", ylabel="Frequency")    
mean(coarse_est_vec)
std(coarse_est_vec)

histogram(fine_est_vec, title="Fine estimation of l_2", label="l_2 fine estimates", xlabel="l_2 estimate", ylabel="Frequency")    
mean(fine_est_vec)
std(fine_est_vec)
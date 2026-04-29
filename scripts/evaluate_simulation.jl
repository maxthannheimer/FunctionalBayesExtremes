using DrWatson
using Distributions
#@quickactivate :FunctionalBayesExtremes
include(srcdir("FunctionalBayesExtremes_MODULE.jl"))
using .FunctionalBayesExtremes_MODULE
using PrettyTables
""" RSME, BIAS, empirical coverage and interval width calculation for the conditional simulation and the approximation method"""

####################
#set name of the simulation results folder
date_string="2026_04_28"
####################
N_burn_in=1000
quantile_val=0.1
true_param=Parameter(α=0.5, β=1.5, c=3.0)
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
    est_approx_mean=Dict( "β" => [NaN for i in 1:total_simulation_number], 
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_approx_median=Dict( "β" => [NaN for i in 1:total_simulation_number],
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_approx_lower_quantile=Dict( "β" => [NaN for i in 1:total_simulation_number],
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_approx_upper_quantile=Dict( "β" => [NaN for i in 1:total_simulation_number],
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number]) 
    est_approx_interval_width=Dict( "β" => [NaN for i in 1:total_simulation_number],
                        "c" => [NaN for i in 1:total_simulation_number], 
                        "α" => [NaN for i in 1:total_simulation_number])
    est_approx_empirical_coverage=Dict( "β" => [NaN for i in 1:total_simulation_number], 
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
    RMSE_est_approx_mean=Dict( "β" => NaN, 
                        "c" => NaN, 
                        "α" => NaN)
    RMSE_est_approx_median=Dict( "β" => NaN,
                        "c" => NaN, 
                        "α" => NaN)
    BIAS_approx_mean=Dict( "β" => NaN, 
                        "c" => NaN, 
                        "α" => NaN)
    BIAS_approx_median=Dict( "β" => NaN, 
                        "c" => NaN, 
                        "α" => NaN)
    approx_interval_width=Dict( "β" => NaN,
                        "c" => NaN, 
                        "α" => NaN)
    approx_empirical_coverage=Dict( "β" => NaN,
                        "c" => NaN, 
                        "α" => NaN) 


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
    for key in keys(param_res_dict)          
                    est_cond_sim_mean[key][file_number]=mean(param_res_dict[key][N_burn_in:end])
                    est_cond_sim_sd[key][file_number]=std(param_res_dict[key][N_burn_in:end])
                    est_cond_sim_median[key][file_number]=median(param_res_dict[key][N_burn_in:end])
                    est_cond_sim_lower_quantile[key][file_number]=quantile(param_res_dict[key][N_burn_in:end], quantile_val)
                    est_cond_sim_upper_quantile[key][file_number]=quantile(param_res_dict[key][N_burn_in:end], 1-quantile_val)
                    est_approx_mean[key][file_number]=mean(param_res_dict_approx[key][N_burn_in:end])
                    est_approx_median[key][file_number]=median(param_res_dict_approx[key][N_burn_in:end])
                    est_approx_lower_quantile[key][file_number]=quantile(param_res_dict_approx[key][N_burn_in:end], quantile_val)
                    est_approx_upper_quantile[key][file_number]=quantile(param_res_dict_approx[key][N_burn_in:end], 1-quantile_val)

                    est_cond_sim_empirical_coverage[key][file_number]= ((est_cond_sim_lower_quantile[key][file_number] .<= true_param_dict[key] ) .* (est_cond_sim_upper_quantile[key][file_number] .>= true_param_dict[key] ))*1.0
                    #est_cond_sim_empirical_coverage[key][number_parallel_int,number_sim_int]=sum( (est_cond_sim_lower_quantile[key][number_parallel_int,number_sim_int] .<= true_param[key] ) .* (est_cond_sim_upper_quantile[key][number_parallel_int,number_sim_int] .>= true_param[key] ) )*1.0
                  
                    (lower_quantile, upper_quantile)=quantile.(Normal(est_cond_sim_mean[key][file_number],est_cond_sim_sd[key][file_number]), [quantile_val, 1-quantile_val])
                    est_cond_sim_normal_coverage[key][file_number]=( lower_quantile .<= true_param_dict[key] ) .* (upper_quantile .>= true_param_dict[key] )*1.0
                    est_approx_empirical_coverage[key][file_number]=( (est_approx_lower_quantile[key][file_number] .<= true_param_dict[key] ) .* (est_approx_upper_quantile[key][file_number] .>= true_param_dict[key] ) )*1.0
                    #cond_sim_empirical_coverage[key][number_parallel_int,number_sim_int]=sum( (est_cond_sim_lower_quantile[key][number_parallel_int,number_sim_int] .<= true_param[key] ) && (est_cond_sim_upper_quantile[key][number_parallel_int,number_sim_int] .>= true_param[key] ) )*1.0
                    #approx_empirical_coverage[key][number_parallel_int,number_sim_int]=sum( (est_approx_lower_quantile[key][number_parallel_int,number_sim_int] .<= true_param[key] ) && (est_approx_upper_quantile[key][number_parallel_int,number_sim_int] .>= true_param[key] ) )*1.0
                    est_cond_sim_interval_width[key][file_number]=est_cond_sim_upper_quantile[key][file_number]-est_cond_sim_lower_quantile[key][file_number]
                    est_approx_interval_width[key][file_number]=est_approx_upper_quantile[key][file_number]-est_approx_lower_quantile[key][file_number]
    end

    for key in keys(est_approx_mean)
                RMSE_cond_sim_mean[key]=round(mean((est_cond_sim_mean[key].-true_param_dict[key]).^2)^0.5 , digits=3)
                RMSE_cond_sim_median[key]=round(mean((est_cond_sim_median[key].-true_param_dict[key]).^2)^0.5 , digits=3)
                BIAS_cond_sim_mean[key]=round(mean(est_cond_sim_mean[key].-true_param_dict[key]) , digits=3)
                BIAS_cond_sim_median[key]=round(mean(est_cond_sim_median[key].-true_param_dict[key]), digits=3)

                cond_sim_empirical_coverage[key]=round(mean(est_cond_sim_empirical_coverage[key]), digits=3)
                cond_sim_normal_coverage[key]=round(mean(est_cond_sim_normal_coverage[key]), digits=3)
                cond_sim_interval_width[key]=round(mean((est_cond_sim_interval_width[key]).^2)^0.5 , digits=3)
            
                RMSE_est_approx_mean[key]=round(mean((est_approx_mean[key].-true_param_dict[key]).^2)^0.5 , digits=3)
                RMSE_est_approx_median[key]=round(mean((est_approx_median[key].-true_param_dict[key]).^2)^0.5 , digits=3)
                BIAS_approx_mean[key]=round(mean(est_approx_mean[key].-true_param_dict[key]) , digits=3)
                BIAS_approx_median[key]=round(mean(est_approx_median[key].-true_param_dict[key]) , digits=3)

                approx_empirical_coverage[key]=round(mean(est_approx_empirical_coverage[key]), digits=3)
                approx_interval_width[key]=round(mean((est_approx_interval_width[key]).^2)^0.5 , digits=3)
    end
end

    #println("RMSE_est_approx_median and median")
    RMSE_est_approx_mean
    RMSE_est_approx_median
    #println("BIAS_approx_mean and median")
    BIAS_approx_mean
    BIAS_approx_median
    #println("approx_empirical_coverage and interval width")
    approx_empirical_coverage
    approx_interval_width
    #println("RMSE_cond_sim_mean and median")
    RMSE_cond_sim_mean
    RMSE_cond_sim_median
    #println("BIAS_cond_sim_mean and median")
    BIAS_cond_sim_mean
    BIAS_cond_sim_median
    #println("cond_sim_empirical_coverage, normal coverage and interval width")
    cond_sim_empirical_coverage
    cond_sim_normal_coverage
    cond_sim_interval_width



   
""" print results in a pretty table"""   
key_array=["β", "c", "α"]
pretty_table_data=Dict{String,Vector{Any}}()
    for key in key_array
        pretty_table_data[key]=[true_param_dict[key],RMSE_est_approx_mean[key], RMSE_est_approx_median[key], BIAS_approx_mean[key], BIAS_approx_median[key], approx_empirical_coverage[key], 
                                RMSE_cond_sim_mean[key], RMSE_cond_sim_median[key], BIAS_cond_sim_mean[key], BIAS_cond_sim_median[key], cond_sim_empirical_coverage[key], cond_sim_normal_coverage[key]]
    end
data=hcat([pretty_table_data[key] for key in key_array]...)'
column_labels =[["Parameter", "RMSE_Mean", "RMSE_Median", "BIAS_Mean", "BIAS_Median", "Emp_Coverage", "RMSE_Mean", "RMSE_Median", "BIAS_Mean", "BIAS_Median", "Emp_Coverage", "Normal_Coverage"],
                [ "", "approx", "", "", "", "", "cond", "", "", "", "", ""]]
pretty_table(data; column_labels)




"""Plotting α,β and c for different realisations of the MCMC and MCMC_approx chains"""

using Plots
file_number=0



    file_number=file_number+1
    file_number=48
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
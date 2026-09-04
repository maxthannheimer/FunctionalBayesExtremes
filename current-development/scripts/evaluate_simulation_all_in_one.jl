
include(joinpath(@__DIR__,"..", "src", "FunctionalBayesExtremes_MODULE.jl"))
using .FunctionalBayesExtremes_MODULE
using PrettyTables
using JLD2
using Distributions
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
total_simulation_number = size(readdir(joinpath(@__DIR__, "..", "data", "exp_raw", date_string)), 1)
#total_simulation_number=size(readdir(datadir("exp_raw",date_string)),1)


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

#data_dict_tmp=load(joinpath(@__DIR__, "..", "data", "exp_raw", date_string, file_string[1]))
for sim_data_name in sim_data_names                    

#here starts the actual calculation
for file_number in 1:size(readdir(joinpath(@__DIR__, "..", "data", "exp_raw", date_string)),1)
#file_number=1 
file_string=readdir(joinpath(@__DIR__, "..", "data", "exp_raw", date_string))
    data_dict_tmp=load(joinpath(@__DIR__, "..", "data", "exp_raw", date_string, file_string[file_number]))
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










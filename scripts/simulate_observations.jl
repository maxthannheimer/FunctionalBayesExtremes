using DrWatson
@quickactivate :FunctionalBayesExtremes
using Distributions

#@quickactivate "FunctionalBayesExtremes"
#using .FunctionalBayesExtremes
""" import function for simulating observations of the Pareto process from Pareto_simu.jl """
#simulate_pareto_process()=FunctionalBayesExtremes.simulate_pareto_process()


""" Generate some observation data and save it to the data directory. """

""" define parameters """
allparams = Dict(
    "α" => [0.5, 1.5],
    "β" => [0.5, 1.5],
    "c" => [0.5, 1.5],
)

""" create parameter dicts for all combinations of parameters """
dicts=dict_list(allparams)


grid=default_Grid()
num_sim=10 #-> 50000
num_runs=200
# i,d= first(enumerate(dicts))
# param=Parameter(α=d["α"], β=d["β"], c=d["c"])



for (i,d) in enumerate(dicts)
    param=Parameter(α=d["α"], β=d["β"], c=d["c"])
    println(param)
    @time( (sim_data, obs_data, obs_x0) = FunctionalBayesExtremes.simulate_pareto_process(param=param, grid=grid, num_runs=num_runs, num_sim=num_sim))
    sim_res = Dict("param" => param, "sim_data" => sim_data, "obs_data" => obs_data, "obs_x0" => obs_x0)
    safesave(datadir("observations", "obs_$(num_sim)_$(i).jld2"), sim_res)
end


#using Distributions
quantile_lvl=0.95
#threshold=10.0
x_1_row=23
x_2_row=47




for (i,d) in (enumerate(dicts))


x_1_coord=grid.coord_fine[x_1_row,:]
x_2_coord=grid.coord_fine[x_2_row,:]
sim_res=load(datadir("observations", "obs_$(num_sim)_$(i).jld2"))["sim_data"]
param=load(datadir("observations", "obs_$(num_sim)_$(i).jld2"))["param"]
x_1_quantile_thresh = quantile(sim_res[:,x_1_row], quantile_lvl)
x_2_quantile_thresh = quantile(sim_res[:,x_2_row], quantile_lvl)
sum_x_1_and_x_2_exceedances=0
sum_x_1_exceedances=0
sum_x_2_exceedances=0
for i in 1:num_sim
    if sim_res[i,x_1_row]>x_1_quantile_thresh
        sum_x_1_exceedances += 1
        if sim_res[i,x_2_row]>x_2_quantile_thresh 
        sum_x_1_and_x_2_exceedances += 1
        end 
    end
    if sim_res[i,x_2_row]>x_2_quantile_thresh
        sum_x_2_exceedances += 1
    end
end

sum_x_1_and_x_2_exceedances
sum_x_1_exceedances
sum_x_2_exceedances
emp_tail_dep_coeff=sum_x_1_and_x_2_exceedances/sum_x_1_exceedances #empirical tail dependence coefficient


vario_value=FunctionalBayesExtremes.vario(coord=x_2_coord-x_1_coord, param=param)
true_tail_dep_coeff=2*(1-cdf.(Normal(0,1),sqrt(vario_value)/sqrt(2))) #theoretical tail dependence coefficient
println("parameters: ", "α=", param.α, ", β=", param.β, ", c=", param.c)
println("empirical tail dependence coefficient: ", emp_tail_dep_coeff)
println("theoretical tail dependence coefficient: ", true_tail_dep_coeff)

println("exceedances at x1: ", sum_x_1_exceedances)
println("exceedances at x2: ", sum_x_2_exceedances)

#println("quantile at x1: ", x_1_quantile_thresh)
#println("quantile at x2: ", x_2_quantile)
end
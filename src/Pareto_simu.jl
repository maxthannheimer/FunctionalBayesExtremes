""" Simulation of Pareto process realizations using MCMC approach from Dombry et al 2015 
    simulate chain with pareto process as stationary distribution """
function simu_specfcts_MCMC(;param::Parameter, grid::Grid, num_runs::Int)::Vector{Float64}
    #sample num_runs+1 many realizations of exp(1/α[G(s)-G(x0)-γ(s-x0)])
    tmp = r_log_gaussian(param=param, grid=grid, num_sim=num_runs+1) 
    old_value = tmp[1]
    #use samples in MCMC approach to converge to stationary distrib W^(r)
    for trial in 1:num_runs
            proposal = tmp[trial+1]
            acceptance_rate = min(1,mean(proposal)^param.α/mean(old_value)^param.α)   
            if (rand()< acceptance_rate)
                old_value=proposal
            end
    end
    #normalize sample and multiply pareto intesity to get pareto process sample
    old_value = old_value/mean(old_value)
    return old_value*=(1/(1-rand()))^(1/param.α)
end

#TODO calculate extremal dependence/upper tail dependence coefficient

#simu_specfcts_MCMC(param=param, grid=grid, num_runs=1000)

function simulate_pareto_process(;param::Parameter, grid::Grid, num_runs::Int, num_sim::Int)#::Vector{Vector{Float64}}
    sim_data = [simu_specfcts_MCMC(param=param, grid=grid, num_runs=num_runs) for i in 1:num_sim]
    sim_data=reduce(hcat,sim_data)' #just make vector of vectors a matrix (same below for observations)
    observation_data=reduce(hcat,[sim_data[i,grid.rows_coord_coarse] for i in 1:num_sim])' #first argument is number of sim, second is coordinate
    observation_x0=vec(reduce(hcat,[sim_data[i,grid.x0] for i in 1:num_sim]))    
    return sim_data,observation_data,observation_x0
end



#sim_data, obs_data, obs_x0 = simulate_pareto_process(param=param, grid=grid, num_runs=num_runs, num_sim=num_sim)
#obs=Observation(param=param, grid=grid, num_runs=1000, num_sim=100)

#obs.sim_data
#obs.obs_data
#obs.obs_x0
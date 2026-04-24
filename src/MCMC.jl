export MCMC, MCMC_approx


#this function proposes a new value uniformly distributed in a 2 eps window in between min_val and max_val arround the old_param
function uniform_proposal(old_param,eps,min_val,max_val)
    if (old_param>min_val+eps && old_param<max_val-eps)
        new_param=rand(Uniform(old_param-eps,old_param+eps))
        old_interval=2*eps
    elseif (old_param<min_val+eps)
        new_param=rand(Uniform(min_val,old_param+eps))#,old_param+eps-min_val
        old_interval=old_param+eps-min_val
    else
        new_param=rand(Uniform(old_param-eps,max_val)) # max_val-old_param-eps
        old_interval=max_val-old_param-eps
    end
    if (new_param>min_val+eps && new_param<max_val-eps)
        new_interval=2*eps
    elseif (new_param<min_val+eps)
        #new_param=rand(Uniform(min_val,old_param+eps))#,old_param+eps-min_val
        new_interval=new_param+eps-min_val
    else
        # =rand(Uniform(old_param-eps,max_val)) # max_val-old_param-eps
        new_interval=max_val-new_param-eps
    end
    return (new_param,old_interval,new_interval)
end

function gaussian_proposal(old_param,eps)
    return exp(rand(Normal(log(old_param),eps)))
end

function parameter_update(;param_old=param_old,param_new=param_new,log_likelihood_old=log_likelihood_old,log_likelihood_new=log_likelihood_new)
    #calculate MCMC acceptance rate a
    if !isfinite(log_likelihood_new)
        a=0
    elseif log_likelihood_new<log_likelihood_old
        a=exp(log_likelihood_new-log_likelihood_old)
    else
        a=1
    end
    rand_sample=rand()
    #accept new value with probability set as acceptance rate a
    if (rand_sample<a)
        #we update param aka return the new ones
        return(param_new,log_likelihood_new)
    else
        #else we just keep the old values and return them
        return(param_old,log_likelihood_old)
    end
end

function log_likehood_log_gauss_1d_non_normalized(x,mu,sigma)
    (-0.5*(log(x)-mu)^2/sigma^2)
end



function MCMC(;N_MCMC::Int,observation::Observation,threshold::Float64,param::Parameter, grid::Grid,N_est_c::Int,N_cond_sim::Int,N_est_d::Int)
    num_obs=size(observation.obs_x0,1)

    param_vec = [Parameter(α=NaN, β=NaN, c=NaN) for i=1:N_MCMC+1]
    param_vec[1]=param
    number_exceed_vec = [NaN for i=1:N_MCMC]
    log_likelihood_vec = [NaN for i=1:N_MCMC]
    res_ell_X_vec = [[NaN for i in 1:num_obs] for j in 1:N_MCMC] 
    
    for trial in 1:N_MCMC


        (exceedance_observation, res_ell_X_vec[trial]) = exceed_cond_sim(num_runs=N_cond_sim, observation=observation, threshold=threshold, param=param_vec[trial], grid=grid )

        #propose new params
        eps_beta=0.05 # 0.05 -> 0.025
        eps_c=0.1 #0.1 -> 0.05
        eps_alpha=0.1 #0.1 -> 0.05
        beta_eps,old_interval,new_interval=uniform_proposal(param_vec[trial].β,eps_beta,0.0,2.0)
        param_eps=Parameter(c=gaussian_proposal(param_vec[trial].c,eps_c),β=beta_eps,α=gaussian_proposal(param_vec[trial].α,eps_alpha))

      
        l1=l_1_fun(param=param_vec[trial], grid=grid,exceedance_observation=exceedance_observation)   
        l2=l_2_fun(param=param_vec[trial],grid=grid,N_est_c=N_est_c,exceedance_observation=exceedance_observation)
        l3=l_3_fun(exceedance_observation=exceedance_observation, param=param_vec[trial], threshold=1.0)
        l4=l_4_fun(exceedance_observation=exceedance_observation, threshold=1.0, param=param_vec[trial], grid=grid ,N_est_d=N_est_d)





        prior = log_likehood_log_gauss_1d_non_normalized(param_vec[trial].c,0.0,1.5)+log_likehood_log_gauss_1d_non_normalized(param_vec[trial].α,0.0,1.0)
        log_likelihood_old=sum([l1,l2,l3,l4,prior,-log(old_interval)]) #-log interval length of uniform proposal for beta
      
        l1=l_1_fun(param=param_eps, grid=grid,exceedance_observation=exceedance_observation)   
        l2=l_2_fun(param=param_eps,grid=grid,N_est_c=N_est_c,exceedance_observation=exceedance_observation)
        l3=l_3_fun(exceedance_observation=exceedance_observation, param=param_eps, threshold=1.0)
        l4=l_4_fun(exceedance_observation=exceedance_observation, threshold=1.0, param=param_eps, grid=grid ,N_est_d=N_est_d)




        prior = log_likehood_log_gauss_1d_non_normalized(param_eps.c,0.0,1.5)+log_likehood_log_gauss_1d_non_normalized(param_eps.α,0.0,1.0)

        log_likelihood_new=sum([l1,l2,l3,l4,prior,-log(new_interval)]) #-log interval length of uniform proposal for beta

        #MCMC update of param according to acceptance rate calculated with old and new likelihoods
        #update and safe param after MCMC step
        param_vec[trial+1],log_likelihood_vec[trial] =parameter_update(param_old=param_vec[trial],param_new=param_eps,log_likelihood_old=log_likelihood_old,log_likelihood_new=log_likelihood_new)
        number_exceed_vec[trial]=size(exceedance_observation.obs_x0,1)
    end
    Dict( "param" => param_vec, "Number of exceedance"=> number_exceed_vec, "log_likelihood" => log_likelihood_vec, "res_ell_X" => res_ell_X_vec) # "log_likelihood_diff" => log_likelihood_diff_vec, "log_likelihood_diff_short" => log_likelihood_diff_vec_short)
end  



function MCMC_approx(;N_MCMC::Int,observation::Observation,threshold::Float64,param::Parameter, grid::Grid,N_est_c::Int,N_cond_sim::Int,N_est_d::Int)
    num_obs=size(observation.obs_x0,1)

    param_vec = [Parameter(α=NaN, β=NaN, c=NaN) for i=1:N_MCMC+1]
    param_vec[1]=param
    number_exceed_vec = [NaN for i=1:N_MCMC]
    log_likelihood_vec = [NaN for i=1:N_MCMC]
    res_ell_X_vec = [[NaN for i in 1:num_obs] for j in 1:N_MCMC] 
    
    for trial in 1:N_MCMC


        (exceedance_observation, res_ell_X_vec[trial]) = exceed_cond_sim_approx(observation=observation, threshold=threshold)

        #propose new params
        eps_beta=0.05 # 0.05 -> 0.025
        eps_c=0.1 #0.1 -> 0.05
        eps_alpha=0.1 #0.1 -> 0.05
        beta_eps,old_interval,new_interval=uniform_proposal(param_vec[trial].β,eps_beta,0.0,2.0)
        param_eps=Parameter(c=gaussian_proposal(param_vec[trial].c,eps_c),β=beta_eps,α=gaussian_proposal(param_vec[trial].α,eps_alpha))

      
        l1=l_1_fun(param=param_vec[trial], grid=grid,exceedance_observation=exceedance_observation)   
        l2=l_2_fun_approx(param=param_vec[trial],grid=grid,N_est_c=N_est_c,exceedance_observation=exceedance_observation)
        l3=l_3_fun(exceedance_observation=exceedance_observation, param=param_vec[trial], threshold=1.0)




       
       
        prior = log_likehood_log_gauss_1d_non_normalized(param_vec[trial].c,0.0,1.5)+log_likehood_log_gauss_1d_non_normalized(param_vec[trial].α,0.0,1.0)
        log_likelihood_old=sum([l1,l2,l3,prior,-log(old_interval)]) #-log interval length of uniform proposal for beta
      
        l1=l_1_fun(param=param_eps, grid=grid,exceedance_observation=exceedance_observation)   
        l2=l_2_fun_approx(param=param_eps,grid=grid,N_est_c=N_est_c,exceedance_observation=exceedance_observation)
        l3=l_3_fun(exceedance_observation=exceedance_observation, param=param_eps, threshold=1.0)





        prior = log_likehood_log_gauss_1d_non_normalized(param_eps.c,0.0,1.5)+log_likehood_log_gauss_1d_non_normalized(param_eps.α,0.0,1.0)
        log_likelihood_new=sum([l1,l2,l3,prior,-log(new_interval)]) #-log interval length of uniform proposal for beta

        #MCMC update of param according to acceptance rate calculated with old and new likelihoods
        #update and safe param after MCMC step
        param_vec[trial+1],log_likelihood_vec[trial] =parameter_update(param_old=param_vec[trial],param_new=param_eps,log_likelihood_old=log_likelihood_old,log_likelihood_new=log_likelihood_new)
        number_exceed_vec[trial]=size(exceedance_observation.obs_x0,1)
    end
    Dict( "param" => param_vec, "Number of exceedance"=> number_exceed_vec, "log_likelihood" => log_likelihood_vec, "res_ell_X" => res_ell_X_vec) # "log_likelihood_diff" => log_likelihood_diff_vec, "log_likelihood_diff_short" => log_likelihood_diff_vec_short)
end  
library(truncdist)

reader = function(file_path = NULL, file_prefix = 'Turbine', file_suffix = '_2017.csv', n_turbines = 66){
  data_frames = vector("list", n_turbines)
  for (i in 1:n_turbines){
    data_frames[[i]] = read.csv(paste0(file_path, file_prefix, i, file_suffix))  
  }
  return(data_frames)
}

thin_data = function(df, start_idx=1, thinning_number=15){
  if (start_idx > thinning_number){
    stop("starting index should not be greater than the thinning number")
  }
  end_idx = min(sapply(df, nrow))
  idx = seq(start_idx, end_idx, thinning_number)
  return(lapply(df, function(x) x[idx,]))
}

computeThinningNumber = function(X){
  thinningNumber = max(apply(X,2,function(col) 
    min(which(c(1,abs(stats::pacf(col, plot = FALSE)$acf[,1,1])) <= (2/sqrt(nrow(X)))))))
  return(thinningNumber)
}

logistic_powercurve = function(wind_speed, theta, temperature, eta){
  if (is.null(dim(theta))){
    return(1/(1 + exp(-(theta[1] + (eta[1]*temperature))*(wind_speed - (theta[2] + (eta[2]*temperature))))))
  } else {
    if (ncol(wind_speed) != nrow(theta)){
      stop("dimension mismatch between theta and wind_speed")
    }
    return(t(1/(1 + exp(-(theta[,1] + (eta[,1]*t(temperature)))*(t(wind_speed) - (theta[,2] + (eta[,2]*t(temperature))))))))
  }
}

# Norm_vector
norm_vec = function( x , y = rep(0,length(x)) ) {
  norm_vector = sqrt(sum((x-y)^2))
  
  if (norm_vector == Inf){
    norm_vector = 1e+308
  } else {
    norm_vector = norm_vector
  } 
  return(norm_vector)
}

likelihood_ratio_theta = function(power, wind_speed, temperature, eta, theta_new, theta_old, sigma_sq) {
  pred_power_new = logistic_powercurve(wind_speed, theta_new, temperature, eta)
  pred_power_old = logistic_powercurve(wind_speed, theta_old, temperature, eta)
  loglik_diff = (sum(pred_power_old**2) - sum(pred_power_new**2) + (2*sum(power*(pred_power_new - pred_power_old))))/(2*sigma_sq)
  #loglik_diff = sum(pnorm(power, mean = pred_power_new, sd = sqrt(sigma_sq), log.p = TRUE)) - sum(pnorm(power, mean = pred_power_old, sd = sqrt(sigma_sq), log.p = TRUE))
  return(loglik_diff)
}

likelihood_ratio_eta = function(power, wind_speed, temperature, theta, eta_new, eta_old, sigma_sq) {
  pred_power_new = logistic_powercurve(wind_speed, theta, temperature, eta_new)
  pred_power_old = logistic_powercurve(wind_speed, theta, temperature, eta_old)
  loglik_diff = (sum(pred_power_old**2) - sum(pred_power_new**2) + (2*sum(power*(pred_power_new - pred_power_old))))/(2*sigma_sq)
  #loglik_diff = sum(pnorm(power, mean = pred_power_new, sd = sqrt(sigma_sq), log.p = TRUE)) - sum(pnorm(power, mean = pred_power_old, sd = sqrt(sigma_sq), log.p = TRUE))
  return(loglik_diff)
}

# Sampling from multivariate Guassian distribution
rmvt.Rue.sample = function(Q, b){
  
  # Goal:
  # Sample from beta = N[mu, Sigma]
  # such that
  # mu = solve(Q)%*%b # p dim vector
  # Sigma = solve(Q) # p by p matrix
  
  # Q : p by p matrix, Precision matrix
  # b : p dim vector
  
  # Useful 1. when n >>> p (i.e. Number of data is more than covariate)
  #        2. Possibly, We want to utilize precision structure
  #        3. Use cholesky and avoiding inverting precision matrix by solving three linear equations
  p = dim(Q)[1]
  
  # Step 1
  L = t(chol(Q))
  
  # Step 2
  z = rnorm(n = p ,mean = 0,sd = 1)
  
  # Step 3
  y = solve(t(L), z)
  
  # Step 4
  v = solve(L, b)
  
  # Step 5
  theta = solve(t(L), v)
  
  # Step 6
  beta = y + theta
  
  return(beta)
  
}


# Ingredient of slice sampler
updated.x = function(old.x, a,b){
  
  # library(truncdist)
  # Goal: sampling from truncated inverse gamma distribution
  # x ~ pi(x) \propto IG[x|a,b] * g(x) such that g(x) = x/(1+x)
  
  # Step 1 : u|x ~ unif(0,g(x))
  u = runif(n = 1, min = 0, max = old.x/(1+old.x)) 
  
  # Step 2 : x|u ~ Trucated inverse gamma distribution
  x = 1/rtrunc(n=1, spec = "gamma", shape = a, rate = b, a = 0, b = (1-u)/u)
  
  return(x)
}

bhm_mcmc_computation = function(power, wind_speed, temperature, terrain, rated_power=100, burn=10000, nmc=10000, thin=30) {
  
  n_turbines = nrow(terrain)
  n_terrain = ncol(terrain)
  
  terrain_scaler = list('location'=rep(0, n_terrain), 'scale'=rep(1, n_terrain))
  terrain_scaler$location = apply(terrain,2,mean)
  terrain_scaler$scale = apply(terrain,2,sd)
  terrain = sapply(c(1:n_terrain), function(j) (terrain[,j] - terrain_scaler$location[j])/terrain_scaler$scale[j])  
  
  
  normalized_power = power/rated_power
  
  temperature_scaler = list('location'=mean(temperature), 'scale'=sd(temperature))
  temperature = (temperature - temperature_scaler$location)/temperature_scaler$scale
  
  n_obs = nrow(wind_speed) 
  
  n_pc_params = 2 #Number of parameters in the pwoer curve
  
  # Bayesian Hierarchical Model
  # Setting
  S = burn + nmc
  
  # Make rooms for parameters
  # Model parameters & measurement error
  theta = array(0, dim = c(n_turbines, n_pc_params, S)) #A tensor to store theta samples
  alpha = matrix(0, nrow = n_pc_params, ncol = S) #A matrix to store alpha samples
  beta = array(0,dim = c(n_terrain, n_pc_params, S)) #A tensor to store beta samples
  sigma_theta_sq = matrix(0, nrow = n_pc_params, ncol = S) #A matrix to store noise variance samples of theta
  tau = matrix(0, nrow = n_pc_params, ncol = S) #A matrix to store tau samples
  gam = matrix(0, nrow = n_pc_params, ncol = S) #A matrix to store gamma samples
  eta = array(0, dim = c(n_turbines, n_pc_params, S)) #A tensor to store eta samples
  sigma_eta_sq = matrix(0, nrow = n_pc_params, ncol = S) #A matrix to store noise variance samples of eta
  sigma_epsilon_sq = rep(0,S) #A vector to store noise variance samples
  
  
  # Decide initial values
  #theta_init = c(1.1, -0.01,8, 0.005)
  theta_init = c(1, 8)
  theta[,,1] = matrix(theta_init, nrow = n_turbines, ncol = n_pc_params, byrow = T)
  alpha[,1] = theta_init
  beta[,,1] = 1
  sigma_theta_sq[,1] = 1
  tau[,1] = 1
  eta[,,1] = matrix(theta_init/100,  nrow = n_turbines, ncol = n_pc_params, byrow = T)
  gam[,1] = theta_init/100
  sigma_eta_sq[,1] = 1
  #eta_proposal_width = c(0.5, 1.0)
  sigma_epsilon_sq[1] = 1
  
  # One vector
  one_vec_n_turb = rep(1,n_turbines) ; one_vec_n_terrain = rep(1,n_terrain)
  # Identity matrix
  I_n_turb = diag(one_vec_n_turb) ; I_n_terrain = diag(one_vec_n_terrain)
  
  # Gibbs sampler
  for (s in 1:(S-1)){
    
    #Step 1: Updating theta's sequentially. Sample theta from N-dimensional distribution
    {
      #copy the previous theta sample (s) as the new sample (s+1) beofre sampling theta
      theta_new = theta[,,s]
      for (i in 1:n_pc_params){
        
        #update theta_old for each theta paramater
        theta_old = theta_new
        
        # a. choose an ellipse
        mu_vec = alpha[i,s] + terrain%*%beta[,i,s]
        nu_vec = c(mvtnorm::rmvnorm(n = 1, mean = mu_vec, sigma = sigma_theta_sq[i,s]*I_n_turb))
        
        # b. define a criterion function
        # defined above
        # c. choose a threshold and fix
        log_u = log(runif(n = 1, min = 0, max = 1))
        
        # d. draw an initial proposal
        phi = runif(n = 1, min = -pi, max = pi)
        theta_new[,i] = (theta_old[,i] - mu_vec)*cos(phi) + (nu_vec - mu_vec)*sin(phi) + mu_vec
        
        ll_theta = likelihood_ratio_theta(power=normalized_power, wind_speed=wind_speed,
                                          temperature=temperature, eta=eta[,,s],
                                          theta_new=theta_new, theta_old=theta_old,
                                          sigma_sq=sigma_epsilon_sq[s])
        
        #e.
        if (log_u >= ll_theta){
          
          # define a bracket
          phi_min = -pi ; phi_max = pi
          while (log_u >= ll_theta){
            # shrink the bracket and try a new point
            if (phi > 0){ phi_max = phi } else { phi_min = phi }
            phi = runif(n = 1, min = phi_min, max = phi_max)
            theta_new[,i] = (theta_old[,i] - mu_vec)*cos(phi) + (nu_vec - mu_vec)*sin(phi) + mu_vec
            ll_theta = likelihood_ratio_theta(power=normalized_power, wind_speed=wind_speed,
                                              temperature=temperature, eta=eta[,,s],
                                              theta_new=theta_new, theta_old=theta_old,
                                              sigma_sq=sigma_epsilon_sq[s])
            
          }
          
        }
      }
      
      theta[,,(s+1)]= theta_new
      
    }
    
    # Step 2 : updating alpha
    for (i in 1:n_pc_params){
      term.1 = (1/n_turbines)*t(one_vec_n_turb)%*%(theta[,i,(s+1)] - terrain%*%beta[,i,s])
      term.2 = sqrt(sigma_theta_sq[i,s]/n_turbines)
      alpha[i,s+1] = rnorm(n = 1, mean = term.1, sd = term.2)
    }
    
    # Step 3 : updating beta
    for (i in 1:n_pc_params){
      Q = (1/sigma_theta_sq[i,s])*(t(terrain)%*%terrain + (1/tau[i,s]^2)*diag(n_terrain))
      b = (1/sigma_theta_sq[i,s])*(t(terrain)%*%(theta[,i,(s+1)] -one_vec_n_turb*alpha[i,s+1]))
      beta[,i,(s+1)] = rmvt.Rue.sample(Q = Q, b = b)
    }
    
    
    # Step 4 : updating sigma_theta_sq
    for (i in 1:n_pc_params){
      term.1 = (n_turbines+n_terrain)/2
      term.2 = norm_vec(x = theta[,i,(s+1)], y = one_vec_n_turb*alpha[i,s+1] + terrain%*%beta[,i,(s+1)])^2
      term.3 = t(beta[,i,(s+1)])%*%beta[,i,(s+1)]/(tau[i,s]^2)
      term.4 = (1/2)*(term.2 + term.3)
      sigma_theta_sq[i,s+1] = 1/rgamma(n = 1, shape = term.1, rate = term.4)
    }
    
    
    # Step 5 : updating tau
    for (i in 1:n_pc_params){
      # 6-a : transformation
      omega = tau[i,s]^2
      # 6-b : slice sampler
      a = (n_terrain+1)/2
      b = t(beta[,i,(s+1)])%*%beta[,i,(s+1)]/(2*sigma_theta_sq[i,(s+1)])
      updated.omega = updated.x(old.x = omega,a = a, b = b)
      # 6-c : transformation back
      tau[i,s+1] = sqrt(updated.omega)
    }
    
    # Step 6 : updating eta
    {
      #copy the previous eta sample (s) as the new sample (s+1) before sampling eta
      eta_new = eta[,,s]
      for (i in 1:n_pc_params){
        
        #update eta_old for each eta paramater
        eta_old = eta_new
        
        # a. choose an ellipse
        mu_vec = gam[i,s]*one_vec_n_turb
        nu_vec = c(mvtnorm::rmvnorm(n = 1, mean = mu_vec, sigma = sigma_eta_sq[i,s]*I_n_turb))
        
        # b. define a criterion function
        # defined above
        # c. choose a threshold and fix
        log_u = log(runif(n = 1, min = 0, max = 1))
        
        # d. draw an initial proposal
        phi = runif(n = 1, min = -pi, max = pi)
        eta_new[,i] = ((eta_old[,i] - mu_vec)*cos(phi)) + ((nu_vec - mu_vec)*sin(phi)) + mu_vec
        
        #e.
        ll_eta = likelihood_ratio_eta(power=normalized_power, wind_speed=wind_speed,
                                      temperature=temperature, theta=theta[,,(s+1)],
                                      eta_new=eta_new, eta_old=eta_old,
                                      sigma_sq=sigma_epsilon_sq[s])
        
        if (log_u >= ll_eta){
          # define a bracket
          phi_min = -pi ; phi_max = pi
          while (log_u >= ll_eta){
            # shrink the bracket and try a new point
            if (phi > 0){ phi_max = phi } else { phi_min = phi }
            phi = runif(n = 1, min = phi_min, max = phi_max)
            eta_new[,i] = ((eta_old[,i] - mu_vec)*cos(phi)) + ((nu_vec - mu_vec)*sin(phi)) + mu_vec
            ll_eta = likelihood_ratio_eta(power=normalized_power, wind_speed=wind_speed,
                                          temperature=temperature, theta=theta[,,(s+1)],
                                          eta_new=eta_new, eta_old=eta_old,
                                          sigma_sq=sigma_epsilon_sq[s])
          }
        }
      }
      
      eta[,,(s+1)]= eta_new
      
    }
    
    # Step 7 : updating gamma
    for (i in 1:n_pc_params){
      term.1 = (1/n_turbines)*t(one_vec_n_turb)%*%(eta[,i,(s+1)])
      term.2 = sqrt(sigma_eta_sq[i,s]/n_turbines)
      gam[i,s+1] = rnorm(n = 1, mean = term.1, sd = term.2)
    }
    
    # Step 8 : updating sigma_eta_sq
    for (i in 1:n_pc_params){
      term.1 = (n_turbines-1)/2
      term.2 = sum((eta[,i,(s+1)]- gam[i,(s+1)])**2)
      term.3 = (1/2)*(term.2)
      sigma_eta_sq[i,s+1] = 1/rgamma(n = 1, shape = term.1, rate = term.3)
    }
    
    
    
    # Step 9 : updating sigma_epsilon_sq
    {
      term.1 = n_turbines*n_obs/2
      term.2 = 0.5*sum((normalized_power - logistic_powercurve(wind_speed,
                                                               theta[,,(s+1)], temperature, eta[,,(s+1)]))**2)
      sigma_epsilon_sq[s+1] = 1/rgamma(n = 1, shape = term.1, rate = term.2)
    }
    
    print(s)
  }

  
  output = list(samples=list("theta"=theta,
                             "alpha"=alpha,
                             "beta"=beta,
                             "sigma_theta_sq"=sigma_theta_sq,
                             "tau"=tau,
                             "gam"=gam,
                             "eta"=eta,
                             "sigma_eta_sq"=sigma_eta_sq,
                             "sigma_epsion_sq"=sigma_epsilon_sq),
                data=list("wind_speed"=wind_speed,
                          "scaled_temperature"=temperature,
                          "normalized_power"=normalized_power,
                          "scaled_terrain"=terrain),
                scalers=list("temperature_scaler"=temperature_scaler,
                             "terrain_scaler"=terrain_scaler,
                             "rated_power"=rated_power),
                mcmc=list("burn"=burn, "nmc"=nmc, "thin"=thin))
  
  class(output) = "bhm"
  return(output)
}

predict.bhm = function(bhm_model, test_wind_speed, test_temperature, test_terrain=NULL, ...){
  n_pc_params = dim(bhm_model$samples$sigma_theta_sq)[1]
  if (nrow(test_wind_speed) != nrow(test_temperature)){
    stop("Number of datapoints in wind speed and temperature not the same.")
  }
  
  if(ncol(test_wind_speed) != ncol(test_temperature)){
    stop("ncol wind speed does not match ncol temperature")
  }
  if (!is.null(test_terrain)){
    if (ncol(test_terrain) != ncol(bhm_model$data$scaled_terrain)){
      stop("Number of terrain variables for predictions must be the same as training.")
    }
    if(ncol(test_wind_speed) != nrow(test_terrain)){
      stop("ncol wind speed does not match nrow terrain")
    }  
    for (j in 1:ncol(test_terrain)){
      test_terrain[,j] = (test_terrain[,j] - bhm_model$scalers$terrain_scaler$location[j])/bhm_model$scalers$terrain_scaler$scale[j]
    }
  }
  
  test_temperature = (test_temperature - bhm_model$scalers$temperature_scaler$location)/bhm_model$scalers$temperature_scaler$scale

  mc_index = seq(from = bhm_model$mcmc$burn + 1, to = bhm_model$mcmc$burn + bhm_model$mcmc$nmc, by = bhm_model$mcmc$thin)
  test_pred = matrix(nrow=nrow(test_wind_speed), ncol=ncol(test_wind_speed))
  for (i in 1:ncol(test_wind_speed)){
    test_theta = array(0, dim=c(length(mc_index), n_pc_params))
    test_eta = array(0, dim=c(length(mc_index), n_pc_params))
    if (!is.null(test_terrain)){
      test_theta[,1] = bhm_model$samples$alpha[1,mc_index] + t(test_terrain[i,,drop=FALSE]%*%bhm_model$samples$beta[,1,mc_index])
      test_theta[,2] = bhm_model$samples$alpha[2,mc_index] + t(test_terrain[i,,drop=FALSE]%*%bhm_model$samples$beta[,2,mc_index])
    } else {
      test_theta[,1] = bhm_model$samples$alpha[1,mc_index] 
      test_theta[,2] = bhm_model$samples$alpha[2,mc_index]
    }
    test_eta[,1] = bhm_model$samples$gam[1, mc_index] 
    test_eta[,2] = bhm_model$samples$gam[2, mc_index] 
    test_pred[,i] = rowMeans(logistic_powercurve(
                          wind_speed=matrix(rep(test_wind_speed[,i],length(mc_index)), ncol=length(mc_index)), 
                          theta=test_theta,
                          temperature=matrix(rep(test_temperature[,i], length(mc_index)), ncol=length(mc_index)), 
                          eta=test_eta)
                          )
  }
  
  return(test_pred*bhm_model$scalers$rated_power)
}

get_mc_index = function(bhm_model){
  mc_index = seq(from = bhm_model$mcmc$burn + 1, to = bhm_model$mcmc$burn + bhm_model$mcmc$nmc, by = bhm_model$mcmc$thin)
  return(mc_index)
}



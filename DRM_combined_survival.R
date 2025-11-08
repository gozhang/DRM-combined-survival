### FUNCTION - DRM-BASED UNBIASED SURVIVAL FUNCTION ESTIMATOR USING RIGHT-CENSORED DATA AND LENGTH-BIASED RIGHT-CENSORED DATA ###
### INPUT: rc.dat = LENGTH-BIASED RIGHT-CENSORED DATA (CAN BE GENERATED FROM 'rc.simulate' FUNCTION) ###
### INPUT: lbrc.dat = LENGTH-BIASED RIGHT-CENSORED DATA (CAN BE GENERATED FROM 'lb.c.simulate' FUNCTION) ###
### INPUT: basis_fun = BASIS FUNCTION IN THE DRM FOR LBRC DATA DENSITY FUNCTION RELATIVE TO RC DATA DENSITY FUNCTION ###
### INPUT: method = OPTIMIZATION METHOD FOR NUMERICALLY FITTING THE DRM: EITHER "eqslv" OR "optimization" ###
##############################################################################################################

rclbrcDRM.estimator = function(rc.dat, lbrc.dat, basis_fun, method = "eqslv") {	
  # method = "optimization" means finding the parameters using optim()
  # method = "eqslv" means finding the parameters using equation solver
  
  # Order the right-censored length-biased data based on the observed failure/censoring times
  # Redefine the observed length-biased failure/censoring times and the failure/censoring indicator function
  lbrc.dat = lbrc.dat[order(lbrc.dat$surv.time),]
  times.lb = lbrc.dat$surv.time; indicator.lb = lbrc.dat$delta
  
  # Order the right-censored data based on the observed failure/censoring times
  # Define only the observed failure times and redefine the failure/censoring indicator function
  rc.dat = rc.dat[order(rc.dat$surv.time),]
  times.inc = rc.dat$surv.time[rc.dat$delta==1]; indicator.inc = rc.dat$delta 

  # Confirm whether extra mass point is present if RC data censoring times surpass largest LBRC time 
  # Extra Mass Logical = 1 if True, = 0 Otherwise 
  extra.mass.logical = 0
  if (any(rc.dat$surv.time[rc.dat$delta==0] > max(lbrc.dat$surv.time))) {
    extra.mass.logical = 1
  }

  # Define the combined cohort mass times of the incident failure times and prevalent failure/censoring times
  # Define the failure/censoring indicator function
  # Define the cohort inclusion variable
  # Bind all the variables from both cohorts into a single variable
  times.comb = c(rc.dat$surv.time[rc.dat$delta==1], lbrc.dat$surv.time); indicator.comb = c(rc.dat$delta[rc.dat$delta==1], lbrc.dat$delta)
  cohort.comb = c(rep(1, length(rc.dat$surv.time[rc.dat$delta==1])), rep(0, length(lbrc.dat$surv.time)))
  fit.dat = cbind(times.comb, indicator.comb, cohort.comb)
  fit.dat = data.frame(fit.dat[order(times.comb),])
  times.comb = fit.dat$times; indicator.comb = fit.dat$indicator; cohort.comb = fit.dat$cohort.comb
  
  # Define the combined cohort of the incident/prevalent failure/censoring times, failure/censoring indicator function and cohort inclusion
  # Bind all the variables from both cohorts into a single variable
  obstimes.comb = c(rc.dat$surv.time, lbrc.dat$surv.time); obsindicator.comb = c(rc.dat$delta, lbrc.dat$delta)
  obscohort.comb = c(rep(1, length(rc.dat$surv.time)), rep(0, length(lbrc.dat$surv.time)))
  obs.dat = cbind(obstimes.comb, obsindicator.comb, obscohort.comb)
  obs.dat = data.frame(obs.dat[order(obstimes.comb),])
  obstimes.comb = obs.dat$obstimes.comb; obsindicator.comb = obs.dat$obsindicator.comb; obscohort.comb = obs.dat$obscohort.comb
  
  if (extra.mass.logical == 0) {
    times.comb = sort(unique(times.comb))
  }
  if (extra.mass.logical == 1) {
    times.comb = c(sort(unique(times.comb)), 1e5)
  }
  
  ###################################
  # store the common values to be used later for the sake of computational cost 
  d <- length(basis_fun(x = times.comb[1])) # dim of basis function vector
  n <- length(times.comb)
  basis_fun_times.comb <- matrix(basis_fun(times.comb), nrow = n, ncol = d) # n by d
  ###################################

  theta = matrix(0, nrow = 1, ncol = d)                     # set initial value for DRM parameter
  
  prob.vals.comb = rep(1/(length(times.comb)), length(times.comb))
  tilt_mat <- theta %*% t(basis_fun_times.comb)
  alpha = -log(sum(exp(tilt_mat)*prob.vals.comb)) # set initial value for DRM parameter
  
  q.prob.vals.comb = exp(alpha + tilt_mat)*prob.vals.comb; # q.prob.vals.comb = q.prob.vals.comb/sum(q.prob.vals.comb)
  
  tau = max(times.comb)
  
  pi.val = sum(q.prob.vals.comb*times.comb)/tau
  
  new.prob.vals.comb = rep(0, length(times.comb))					# set starting variable for updated estimated probability mass values
  # new.prob.vals.lb = rep(0, length(prob.vals.lb))
  
  ws = rep(0, length(times.comb))							# set starting variable for w function values
  
  iteration = 0										# set initial iteration value
  eps = 1											# set convergence difference
  conv.crit = 1e-8										# set convergence criterion

  while ((iteration < 2000) & (eps > conv.crit)) {
    
    denoms.p = rep(0, length(obstimes.comb)); denoms.q = rep(0, length(obstimes.comb))
    for (i in 1:length(obstimes.comb)) {
      denoms.p[i] = sum(prob.vals.comb*1*(obstimes.comb[i] <= times.comb))
      denoms.q[i] = sum(q.prob.vals.comb*1*(obstimes.comb[i] <= times.comb))
    }
    
    w1 = rep(0, length(times.comb)) # RC  + LBRC data weights
    w2 = rep(0, length(times.comb)) # LBRC data weights
    
    for (j in 1:length(times.comb)) {
      w1[j] = sum(obs.dat$obscohort.comb*obs.dat$obsindicator.comb*1*(obs.dat$obstimes.comb == times.comb[j])) +
              sum(obs.dat$obscohort.comb*1*(1-obs.dat$obsindicator.comb)*(obs.dat$obstimes.comb <= times.comb[j])*prob.vals.comb[j]/denoms.p) + 
      
              sum((1-obs.dat$obscohort.comb)*obs.dat$obsindicator.comb*1*(obs.dat$obstimes.comb == times.comb[j])) +
              sum((1-obs.dat$obscohort.comb)*(1-obs.dat$obsindicator.comb)*(obs.dat$obstimes.comb <= times.comb[j])*q.prob.vals.comb[j]/denoms.q) +
              ((length(lbrc.dat$surv.time))/pi.val)*(1-(times.comb[j]/tau))*q.prob.vals.comb[j]

      w2[j] = sum((1-obs.dat$obscohort.comb)*obs.dat$obsindicator.comb*1*(obs.dat$obstimes.comb == times.comb[j])) +
              sum((1-obs.dat$obscohort.comb)*(1-obs.dat$obsindicator.comb)*(obs.dat$obstimes.comb <= times.comb[j])*q.prob.vals.comb[j]/denoms.q) +
              ((length(lbrc.dat$surv.time))/pi.val)*(1-(times.comb[j]/tau))*q.prob.vals.comb[j]
    }

    
    if (method == "eqslv"){
      update = find_para_eqslv(times.comb, basis_fun, w1, w2)
    } else if (method == "optimization"){
      update = find_para_opt(times.comb, basis_fun, w1, w2)
    } else {
      stop("Error: method should be either 'eqslv' or 'optimization'.")
    }
    
    alpha = update$alpha; theta = update$theta
    tilt_mat <- theta %*% t(basis_fun_times.comb)

    new.prob.vals.comb = update$p

    new.q.prob.vals.comb = exp(alpha + tilt_mat)*new.prob.vals.comb; # new.q.prob.vals.comb = new.q.prob.vals.comb/sum(new.q.prob.vals.comb)

    pi.val = sum(new.q.prob.vals.comb*times.comb)/tau

    eps = mean(abs(prob.vals.comb - new.prob.vals.comb))
    iteration = iteration + 1

    prob.vals.comb = new.prob.vals.comb
    q.prob.vals.comb = new.q.prob.vals.comb
  }
  
  # Calculate CDF values based on estimated probability mass values
  drops = cumsum(prob.vals.comb)						
  surv = 1 - drops

  surv.est = stepfun(times.comb, c(1, surv))
  return(list(surv.est, alpha, theta))
  
}



##################################################################################################################  
library(nleqslv)
library(rootSolve)
# The following contains the functions that are called in the main function above.
# This code aims to implement the M-step with two approaches: 
# 1) solving the score equation and 2) directly maximizing the likelihood

##################################################################################################################   
# First approach
######################################################################
# The score function for the LBRC+RC DRM
score_LBRC_RC <- function(para, failure_times, basis_fun, w1, w2){ 
  # para is (alpha, theta), where theta has the same dimension as the basis function
  # failure_times is the combined, UNORDERED, failure times from LBRC and RC data, which may contain ties
  # basis_fun is the vector-values basis function
  # w1, w2 are the vectors of weights from the E-step, each having the same length and in the same ORDER as the failure_times
  
  ######################################################################
  # store the common values to be used later for the sake of computational cost 
  d <- length(basis_fun(x = failure_times[1])) # dim of basis function vector
  n <- length(failure_times)
  basis_fun_failure_times <- matrix(basis_fun(failure_times), nrow = n, ncol = d) # n by d
  
  # para = (alpha, theta)
  alpha <- as.vector(para[1])
  theta <- matrix(para[2:(1+d)], ncol = d)
  lambda_1 <- sum(w1-w2)
  lambda_2 <- sum(w2)
  # Note: write in matrix style for easier modification when there are multiple samples
  
  # calculate p  
  tilt_mat <- theta %*% t(basis_fun_failure_times)
  p <- w1/(t(exp(alpha + tilt_mat)) * lambda_2 + lambda_1)
  p <- as.vector(p)
  
  # constraints for (alpha, theta)
  constraint_alpha <- 1 - exp(alpha + tilt_mat) %*% p # an vector of length 1
  q_exp <- (exp(alpha + tilt_mat)) %*% (p*basis_fun_failure_times) # 1 by d
  constraint_theta <- t(basis_fun_failure_times) %*% as.vector(w2) - lambda_2*as.vector(q_exp) # a vector of length 1*d
  constraints <- c(constraint_alpha, constraint_theta)
  # if (min(p) < 0) constraints <- rep(1e8, times = length(constraints))
  
  return(constraints) 
} 

######################################################################
# Solve for para
find_para_eqslv <- function(failure_times, basis_fun, w1, w2){ 
  ######################################################################
  # store the common values to be used later for the sake of computational cost 
  d <- length(basis_fun(x = failure_times[1])) # dim of basis function vector
  n <- length(failure_times)
  basis_fun_failure_times <- matrix(basis_fun(failure_times), nrow = n, ncol = d) # n by d
  lambda_1 <- sum(w1-w2)
  lambda_2 <- sum(w2)
  # Note: write in matrix style for easier modification when there are multiple samples
  
  ###################################
  # Solve the score equations
  obj_fun_marginal <- function(para) return(score_LBRC_RC(para = para, failure_times = failure_times, 
                                                          basis_fun = basis_fun, w1 = w1, w2 = w2))
  # Use the rootsolver() function to solve 
  init_alpha <- rep(0, times = 1)
  init_theta <- rep(0, times = 1*d)
  marginal_opt <- rootsolver(fun = obj_fun_marginal, init_value = c(init_alpha, init_theta))
  
  alpha <- as.vector(marginal_opt$root[1])
  theta <- matrix(marginal_opt$root[(1+1):(1+d*1)], ncol = d)
  # calculate p  
  tilt_mat <- theta %*% t(basis_fun_failure_times)
  p <- w1/(t(exp(alpha + tilt_mat)) * lambda_2 + lambda_1)
  p <- as.vector(p)
  
  is_sum_1 <- abs(1-sum(p)) <= 1e-5
  is_positive <- min(p) >= 0
  
  if ((marginal_opt$root_msg == 0) & is_sum_1 & is_positive){
    # print("Roots for marginal DRM is found!")
  } else {
    print(abs(1-sum(p)))
    print(min(p))
    print(marginal_opt$root_msg)
    stop("Error: roots for marginal DRM are not found!")
  }
  
  return(list(p = p, alpha = alpha, theta = theta)) 
} 

######################################################################
# A self-made optimization function that uses combination of the package rootsolve and nleqslv 
rootsolver <- function(fun, init_value){
  # try multiroot first 
  rootsolver_MR <- multiroot(start = init_value, f = fun, 
                             maxiter = 50, verbose = FALSE, atol = 1e-10)
  
  #if the multiroot succeeds, then no need to try the other root solver
  if ((!is.na(rootsolver_MR$estim.precis)) & (abs(rootsolver_MR$estim.precis) < 1e-6)){
    root <- rootsolver_MR$root
    root_msg <- 0 
  } else { #if the multiroot fails, use the other root solver: nleqslv 
    # try all the combinations of methods and global strategies, and pick the successful one
    test_NL <- testnslv(x = init_value, fn = fun, 
                        method = c("Newton", "Broyden"), 
                        global = c("gline", "pwldog", "hook", "none"),
                        Nrep=0L, title=NULL, control = list(allowSingular = TRUE, maxit = 50))
    test_NL$out$termcd[is.na(test_NL$out$termcd)] <- 4 # specify the NA as a failure#4
    if (1 %in% test_NL$out$termcd){ # 1 in nleqslv means a success
      method_NL <- test_NL$out$Method[test_NL$out$termcd == 1][1]
      global_NL <- test_NL$out$Global[test_NL$out$termcd == 1][1]
      rootsolver_NL <- nleqslv(x = init_value, fn = fun, 
                               global = global_NL, method = method_NL,   
                               control = list(allowSingular = TRUE, maxit = 50, xtol = 1e-12, trace = 0))
      
      root <- rootsolver_NL$x
      root_msg <- rootsolver_NL$termcd - 1 
    } else { # the root can not be found, return the initial values
      root <- init_value
      root_msg <- 10
    }
  }
  
  # return the root and msg (0 means success, 10 means failure)
  return(list(root = root, root_msg = root_msg))
} 


##################################################################################################################   
# Second approach
# The log-likelihood function for the LBRC+RC DRM
log_EL_LBRC_RC <- function(para, failure_times, basis_fun, w1, w2){ 
  # para is (alpha, theta), where theta has the same dimension as the basis function
  # failure_times is the combined, UNORDERED, failure times from LBRC and RC data, which may contain ties
  # basis_fun is the vector-values basis function
  # w1, w2 are the vectors of weights from the E-step, each having the same length and in the same ORDER as the failure_times
  
  ######################################################################
  # store the common values to be used later for the sake of computational cost 
  d <- length(basis_fun(x = failure_times[1])) # dim of basis function vector
  n <- length(failure_times)
  basis_fun_failure_times <- matrix(basis_fun(failure_times), nrow = n, ncol = d) # n by d
  
  # para = (alpha, theta)
  alpha <- as.vector(para[1])
  theta <- matrix(para[2:(1+d)], ncol = d)
  lambda_1 <- sum(w1-w2)
  lambda_2 <- sum(w2)
  # Note: write in matrix style for easier modification when there are multiple samples
  
  # calculate p  
  tilt_mat <- theta %*% t(basis_fun_failure_times)
  p <- w1/(t(exp(alpha + tilt_mat)) * lambda_2 + lambda_1)
  p <- as.vector(p)
  
  # log Empirical Likelihood for (alpha, theta)
  log_lik <- sum(log(p)*w1 + (alpha+t(tilt_mat))*w2)
  
  return(log_lik) 
} 

# Solve for para
find_para_opt <- function(failure_times, basis_fun, w1, w2){ 
  ######################################################################
  # store the common values to be used later for the sake of computational cost 
  d <- length(basis_fun(x = failure_times[1])) # dim of basis function vector
  n <- length(failure_times)
  basis_fun_failure_times <- matrix(basis_fun(failure_times), nrow = n, ncol = d) # n by d
  lambda_1 <- sum(w1-w2)
  lambda_2 <- sum(w2)
  # Note: write in matrix style for easier modification when there are multiple samples
  
  ###################################
  # Maximize the log-EL
  init_alpha <- rep(0, times = 1)
  init_theta <- rep(0, times = 1*d)
  # optim() does minimization by default, so maximizing loglik is the same as minimizing the -loglik
  opt_loglik <- optim(par = c(init_alpha, init_theta), 
                      fn = function(para) -log_EL_LBRC_RC(para=para, failure_times=failure_times, basis_fun=basis_fun, w1=w1, w2=w2),
                      control = list(reltol = 1e-15, maxit = 1000), method = "BFGS")
  # calculate p  
  alpha <- opt_loglik$par[1]
  theta <- matrix(opt_loglik$par[2:(1+d)], ncol = d)
  tilt_mat <- theta %*% t(basis_fun_failure_times)
  p <- w1/(t(exp(alpha + tilt_mat)) * sum(w2) + sum(w1-w2))
  p <- as.vector(p)
  
  is_sum_1 <- abs(1-sum(p)) <= 1e-5
  is_positive <- min(p) >= 0
  
  if ((opt_loglik$convergence == 0) & is_sum_1 & is_positive){
    # print("Roots for marginal DRM is found!")
  } else {
    print(abs(1-sum(p)))
    print(min(p))
    print(opt_loglik$convergence)
    stop("Error: roots for marginal DRM are not found!")
  }
  
  return(list(p = p, alpha = alpha, theta = theta)) 
} 

##################################################################################################################  
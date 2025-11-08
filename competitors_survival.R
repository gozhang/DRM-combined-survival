##################################################################################################################  
# This file contains the functions implementing the competing estimators considered in the paper.

##################################################################################################################  
library(survival)

#' Compute the Kaplan-Meier survival function NPMLE for a set of right-censored failure time data
#'
#' @param fctimes Failure/Censoring times
#' @param fcinds Failure/Censoring indicators
#'
#' return survival step-function estimator

KM.RC.est = function(fctimes, fcinds) {
  surv.obj = survfit(Surv(fctimes, fcinds)~1)
  return(stepfun(surv.obj$time, c(1,surv.obj$surv)))
}

#' Compute the Kaplan-Meier survival function NPMLE for a set of left-truncated right-censored failure time data
#'
#' @param lttimes Left-Truncation times
#' @param fctimes Failure/Censoring times
#' @param fcinds Failure/Censoring indicators
#'
#' return survival step-function estimator

KM.LTRC.est = function(lttimes, fctimes, fcinds) {
  surv.obj = survfit(Surv(lttimes, fctimes, fcinds)~1)
  return(stepfun(surv.obj$time, c(1,surv.obj$surv)))
}


##################################################################################################################  
### FUNCTION - UNBIASED SURVIVAL FUNCTION ESTIMATOR USING LENGTH-BIASED RIGHT-CENSORED DATA ###
### INPUT: lb.dat = LENGTH-BIASED RIGHT-CENSORED DATA (CAN BE GENERATED FROM 'lb.c.simulate' FUNCTION) ###
##################################################################################################

lbrc.emestimator = function(lb.dat) {	
  
  # Reorder LBRC data and rename variables of interest
  lb.dat = lb.dat[order(lb.dat$surv.time),]
  times.lb = lb.dat$surv.time; indicator.lb = lb.dat$delta

  # List time points with positive probability mass
  times.mass = sort(unique(times.lb))
  
  # Set vector for probability masses
  prob.vals = rep(1/(length(times.mass)), length(times.mass))
  
  # Determine the maximum observation for tauhat
  tau = max(times.mass)
  
  # Estimate the probability pi
  pi.val = sum(prob.vals*times.mass)/tau
  
  # Set starting variable for updated estimated probability mass values
  new.prob.vals = rep(0, length(times.mass))					
  
  # Set starting variable for w function values
  ws = rep(0, length(times.mass))							
  
  # Set initial iteration value, convergence difference, convergence criterion
  iteration = 0										
  eps = 1											
  conv.crit = 1e-8						
  
  while ((iteration < 2000) & (eps > conv.crit)) {
    
    denoms = rep(0, length(lb.dat$surv.time))
    
    for (i in 1:length(times.lb)) {denoms[i] = sum(prob.vals*1*(times.lb[i] <= times.mass))}
    
    ws = rep(0, length(times.mass))
    
    for (j in 1:length(times.mass)) {
      ws[j] = sum(lb.dat$delta*1*(lb.dat$surv.time == times.mass[j])) +
        sum(1*(1-lb.dat$delta)*1*(lb.dat$surv.time <= times.mass[j])*prob.vals[j]/denoms) +
        ((length(lb.dat$surv.time))/pi.val)*(1-(times.mass[j]/tau))*prob.vals[j]
      
    }
    
    new.prob.vals = (pi.val/dim(lb.dat)[1])*ws
    
    pi.val = sum(new.prob.vals*times.mass)/tau
    
    eps = mean(abs(prob.vals - new.prob.vals))
    iteration = iteration + 1
    
    prob.vals = new.prob.vals
  }
  
  # Calculate CDF values based on estimated probability mass values
  drops = cumsum(prob.vals)						
  surv = 1 - drops
  
  output = stepfun(times.mass, c(1, surv))
  return(output)
  
}


##################################################################################################################  
### FUNCTION - UNBIASED SURVIVAL FUNCTION ESTIMATOR INCIDENT AND PREVALENT COHORT DATA ###
### INPUT: rc.dat = RIGHT-CENSORED DATA (CAN BE GENERATED FROM 'rc.simulate' FUNCTION) ###
### INPUT: lbrc.dat = PREVALENT COHORT DATA (CAN BE GENERATED FROM 'lb.c.simulate' FUNCTION) ###
######################################################################################

emcomb.estimator = function(rc.dat, lbrc.dat) {	
  
  # Order the right-censored length-biased data based on the observed failure/censoring times
  # Redefine the observed length-biased failure/censoring times and the failure/censoring indicator function
  lbrc.dat = lbrc.dat[order(lbrc.dat$surv.time),]
  times.lb = lbrc.dat$surv.time; indicator.lb = lbrc.dat$delta
  
  # Order the right-censored data based on the observed failure/censoring times
  # Define only the observed failure times and redefine the failure/censoring indicator function
  rc.dat = rc.dat[order(rc.dat$surv.time),]
  times.inc = rc.dat$surv.time[rc.dat$delta==1]; indicator.inc = rc.dat$delta 
  
  # Confirm whether extra mass point is present if incident cohort censoring times surpass largest prevalent cohort time 
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
  
  prob.vals.comb = rep(1/(length(times.comb)), length(times.comb))
  prob.vals.lb = prob.vals.comb[cohort.comb==0]/sum(prob.vals.comb[cohort.comb==0])
  
  tau = max(times.comb)
  
  pi.val = sum(prob.vals.comb*times.comb)/tau
  
  new.prob.vals.comb = rep(0, length(times.comb))					# set starting variable for updated estimated probability mass values
  new.prob.vals.lb = rep(0, length(prob.vals.lb))
  
  ws = rep(0, length(times.comb))							# set starting variable for w function values
  
  iteration = 0										# set initial iteration value
  eps = 1											# set convergence difference
  conv.crit = 1e-8										# set convergence criterion
  
  while ((iteration < 2000) & (eps > conv.crit)) {
    
    denoms = rep(0, length(obstimes.comb))
    for (i in 1:length(obstimes.comb)) {denoms[i] = sum(prob.vals.comb*1*(obstimes.comb[i] <= times.comb))}
    
    ws = rep(0, length(times.comb))
    
    for (j in 1:length(times.comb)) {
      ws[j] = 
        sum(1*(1-obs.dat$obsindicator.comb)*1*(obs.dat$obstimes.comb <= times.comb[j])*prob.vals.comb[j]/denoms) +
        sum(obs.dat$obscohort.comb*obs.dat$obsindicator.comb*1*(obs.dat$obstimes.comb == times.comb[j])) + 
        sum((1-obs.dat$obscohort.comb)*obs.dat$obsindicator.comb*1*(obs.dat$obstimes.comb == times.comb[j])) +
        ((length(lbrc.dat$surv.time))/pi.val)*(1-(times.comb[j]/tau))*prob.vals.comb[j]
    }
    
    new.prob.vals.comb = ws/(dim(rc.dat)[1] + dim(lbrc.dat)[1] + dim(lbrc.dat)[1]*(1-pi.val)/pi.val)
    
    pi.val = sum(new.prob.vals.comb*times.comb)/tau
    
    eps = mean(abs(prob.vals.comb - new.prob.vals.comb))
    iteration = iteration + 1
    
    prob.vals.comb = new.prob.vals.comb
  }
  
  drops = cumsum(prob.vals.comb)						# calculate CDF values based on estimated probability mass values
  surv = 1 - drops
  
  output = stepfun(times.comb, c(1, surv))
  return(output)
  
}

##################################################################################################################  
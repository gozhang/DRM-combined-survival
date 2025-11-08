# R code and demo for the paper  
**Density Ratio Model (DRM) for Multiple Types of Survival Data with Empirical Likelihood**

### Description

This repository provides R code implementing the proposed DRM-based survival function estimation method for combined right-censored (RC) and left-truncated right-censored (LBRC) failure time data.
It also includes implementations of the **competing estimators** (e.g., Kaplan–Meier, NPMLE) considered in the paper, as well as a demo script for reproducing the **simulation studies**.


### Functions

The following are the main R functions for estimating survival functions under the DRM framework and for the competing estimators considered in the paper.


#### **DRM-based estimator**

- **`rclbrcDRM.estimator()`**  
  Estimates the *unbiased survival function* using RC and LBRC data under the **DRM**.  
  Supports two numerical approaches for fitting the DRM:  
  - `"eqslv"` — solving the estimating equations;  
  - `"optimization"` — direct maximization of the empirical likelihood.  

#### **Competing estimators**

Baseline survival estimators used for comparison:

- **`KM.RC.est()`** – Kaplan–Meier NPMLE for RC data.  
- **`KM.LTRC.est()`** – Kaplan–Meier NPMLE for left-truncated RC data.  
- **`lbrc.emestimator()`** – EM-based unbiased survival estimator for LBRC data.  
- **`emcomb.estimator()`** – EM-based unbiased estimator combining RC and LBRC cohorts under the assumption of identical distributions.


### A simulation example

Below is a minimal example showing how to simulate combined RC and LBRC survival data, fit the DRM-based estimator, and compare it with the competing estimators.

```r
# Load required packages
library(survival)
library(nleqslv)
library(rootSolve)

# Function to simulate RC data

### FUNCTION - GENERATE GAMMA DISTRIBUTED RIGHT-CENSORED DATA ###
### INPUTS: n = SAMPLE SIZE, fpar = FAILURE TIME GAMMA c(shp, scl) ###
### cpar = CENSORING TIME EXPONENTIAL LAMBDA ###
################################################
rc.simulate = function(n, fpar, cpar) {     # n = sample size, fpar = fail.time gamma params., cpar = (forward) cens.time exp. param.
  j = 0 
  output = NULL
  while (j < n) {										# iterate until a sample of size n is obtained
    t = rgamma(1, shape=fpar[1], scale=fpar[2])				# generate a failure time 
    cens.t = rexp(1, cpar)    # generate a censoring time
    surv.time = ifelse(t <= cens.t, t, cens.t)    # take the minimum to generate the observed time
    delta = ifelse(t <= cens.t, 1, 0)               # generate indicator value for failure/censoring type
    output = rbind(output, c(surv.time, delta))         # update the output file
    j = j + 1                             # re-iterate
  }
  output = data.frame(output)
  dimnames(output)[[2]] = c("surv.time", "delta")	# give names to output columns for easy access
  return(output)
}

# Function to simulate LBRC data

### FUNCTION - GENERATE GAMMA DISTRIBUTED LENGTH-BIASED RIGHT-CENSORED DATA ###
### INPUTS: n = SAMPLE SIZE, fpar = FAILURE TIME GAMMA c(shp, scl) ###
### cpar = FORWARD CENSORING TIME EXPONENTIAL LAMBDA, R = LOWER BOUND ON UNIFORM ONSET DISTRIBUTION ###
#######################################################################################################
lbrc.simulate = function(n, fpar, cpar, R) {     # n = sample size, fpar = fail.time gamma params., cpar = (forward) cens.time exp. param., R is lower bound on possible onset time from uniform dist.
  j = 0 
  output = NULL
  while (j < n) {										# iterate until a sample of size n is obtained
    onset = runif(1, min = 0, max = R)						# generate an onset time
    t = rgamma(1, shape=fpar[1], scale=fpar[2])				# generate a failure time 
    if (t >= (R-onset)){									# accept the proposed onset/failure time pair if inclusion into the prevalent cohort is satisfied
      ltr = R-onset									# generate the left-truncation time
      obs.t = t										# as inclusion into cohort satisfied, rename failure time 
      cens.t = ltr+rexp(1,cpar)							# generate a forward censoring time (this form of ltr + cens, ensures the censoring time is included in the cohort) 
      surv.time = ifelse(obs.t <= cens.t, obs.t, cens.t)			# take the minimum to generate the observed time
      delta = ifelse(obs.t <= cens.t, 1, 0)					# generate indicator value for failure/censoring type
      output = rbind(output, c(onset, ltr, surv.time, delta))		# update the output file 
      j = j + 1										# re-iterate
    }
  }
  output = data.frame(output)
  dimnames(output)[[2]] = c("onset.time", "ltr.time", "surv.time", "delta")	# give names to output columns for easy access
  return(output)
}

# Simulate combined RC and LBRC data

rc.shape = 0.5
rc.scale = 2
rc.cens.par = 0.2
rc.size = 50
rc.dat = rc.simulate(n = rc.size, fpar = c(rc.shape, rc.scale), cpar = rc.cens.par)

lbrc.shape = rc.shape + 1
lbrc.scale = 2
lbrc.cens.par = 0.07
lbrc.size = 50
lbrc.dat = lbrc.simulate(n = lbrc.size, fpar = c(lbrc.shape, lbrc.scale), cpar = lbrc.cens.par, R = 50)
```

We then run simulation using the simulated data above.

#### 1. DRM approach to estimating survival function and parameters using RC and LBRC data together
```r
# Example basis function used in the paper
basis_fun = function(x) return(log(x))

# Fit the DRM-based estimator by solving the estimating equations
rclbrc.DRM.est_eqslv = rclbrcDRM.estimator(rc.dat, lbrc.dat, basis_fun, method = "eqslv")	

# Or can also directly maximizing the empirical likelihood.  
rclbrc.DRM.est_opt = rclbrcDRM.estimator(rc.dat, lbrc.dat, basis_fun, method = "optimization")	
```

#### 2. Other competitors
```r
# Estimate survival function using RC data only (RC Kaplan-Meier estimate)
rc.est = KM.RC.est(rc.dat$surv.time, rc.dat$delta)
    
# Estimate survival function using LBRC data only (LTRC Kaplan-Meier estimate)
ltrc.est = KM.LTRC.est(lbrc.dat$ltr.time, lbrc.dat$surv.time, lbrc.dat$delta)
    
# Estimate survival function using LBRC data only (LBRC NPMLE)
lbrc.est = lbrc.emestimator(lbrc.dat)
    
# Estimate survival function nonparametrically using RC and LBRC data together (NPMLE, assuming same failure time distributions)
rclbrc.est = emcomb.estimator(rc.dat, lbrc.dat)
```

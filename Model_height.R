model {
  
  #################################################################
  ## 1. LATENT RECOVERY STATES
  #################################################################
  
  ## Recovery state associated with each NDVI observation
  for (o in 1:nObsNDVI) {
    
    recovery.ndvi[o] <-
      1 - exp(-age.ndvi[o] /
                lambda.ndvi[id.ndvi[o]])
    
  }
  
  
  ## Recovery state associated with each LiDAR height observation
  for (h in 1:nObsHeight) {
    
    recovery.height[h] <-
      1 - exp(-age.height[h] /
                lambda.height[id.height[h]])
    
  }
  
  
  #################################################################
  ## 2. NDVI OBSERVATION MODEL
  #################################################################
  
  for (o in 1:nObsNDVI) {
    
    ndvi[o] ~ dnorm(mu.ndvi[o], tau.ndvi)
    
    mu.ndvi[o] <-
      alpha[id.ndvi[o]] +
      gamma[id.ndvi[o]] * recovery.ndvi[o] +
      A[id.ndvi[o]] *
      sin(
        phi +
          ((firemonth.ndvi[o] - 1) * 3.141593 / 6) +
          6.283185 * age.ndvi[o]
      )
    
  }
  
  
  #################################################################
  ## 3. LiDAR VEGETATION HEIGHT OBSERVATION MODEL
  #################################################################
  
  for (h in 1:nObsHeight) {
    
    height[h] ~ dnorm(mu.height[h], tau.height)
    
    ## Hmax = asymptotic vegetation height
    ## H0   = post-fire/baseline height
    
    mu.height[h] <-
      H0[id.height[h]] +
      (Hmax[id.height[h]] -
         H0[id.height[h]]) *
      recovery.height[h]
    
  }
  
  
  #################################################################
  ## 4. GRID-CELL LEVEL PARAMETERS
  #################################################################
  
  for (i in 1:nGrid) {
    
    ## -------------------------------------------------------------
    ## NDVI parameters
    ## -------------------------------------------------------------
    
    alpha[i] ~ dlnorm(alpha.mu, alpha.tau)
    
    gamma[i] ~ dlnorm(gamma.mu[i], gamma.tau)
    
    A[i] ~ dlnorm(A.mu[i], A.tau)
    
    
    ## -------------------------------------------------------------
    ## Height parameters
    ## -------------------------------------------------------------
    
    ## Initial/post-fire vegetation height
    H0[i] ~ dlnorm(H0.mu, H0.tau)
    
    ## Asymptotic vegetation height
    Hmax[i] ~ dlnorm(Hmax.mu[i], Hmax.tau)
    
    
    ## -------------------------------------------------------------
    ## NDVI recovery rate
    ## -------------------------------------------------------------
    
    lambda.ndvi[i] ~
      dlnorm(lambda.ndvi.mu[i],
             lambda.ndvi.tau)
    
    
    ## -------------------------------------------------------------
    ## Hierarchical relationship between recovery rates
    ## -------------------------------------------------------------
    
    ## Height recovery timescale is related to NDVI
    ## recovery timescale:
    ##
    ## log(lambda.height) =
    ##     log(lambda.ndvi) + delta + error
    ##
    ## Therefore:
    ##
    ## lambda.height approximately =
    ##     lambda.ndvi * exp(delta)
    
    log(lambda.height[i]) ~
      dnorm(
        log(lambda.ndvi[i]) + delta.mu,
        delta.tau
      )
    
  }
  
  
  #################################################################
  ## 5. ENVIRONMENTAL REGRESSIONS
  #################################################################
  
  ## Environmental effects on NDVI recovery magnitude
  gamma.mu <- env %*% gamma.beta
  
  ## Environmental effects on NDVI recovery rate
  lambda.ndvi.mu <- env %*% lambda.ndvi.beta
  
  ## Environmental effects on NDVI seasonal amplitude
  A.mu <- env %*% A.beta
  
  ## Environmental effects on asymptotic vegetation height
  Hmax.mu <- env %*% Hmax.beta
  
  
  #################################################################
  ## 6. INTERCEPTS
  #################################################################
  
  ## Mean of log(alpha)
  alpha.mu ~ dnorm(0.15, 10)
  
  ## Mean of log(H0)
  H0.mu ~ dnorm(0, 1)
  
  
  #################################################################
  ## 7. RELATIVE RECOVERY TIMESCALE
  #################################################################
  
  ## Mean log ratio between height and NDVI recovery rates
  ##
  ## delta.mu = 0:
  ##     lambda.height = lambda.ndvi
  ##
  ## delta.mu > 0:
  ##     height recovers more slowly
  ##
  ## delta.mu < 0:
  ##     height recovers more rapidly
  
  delta.mu ~ dnorm(0, 0.25)
  
  delta.tau ~ dgamma(0.01, 0.01)
  
  
  #################################################################
  ## 8. ENVIRONMENTAL REGRESSION COEFFICIENTS
  #################################################################
  
  for (l in 1:nBeta) {
    
    gamma.beta[l] ~ dnorm(0, 0.1)
    
    lambda.ndvi.beta[l] ~ dnorm(0, 0.1)
    
    A.beta[l] ~ dnorm(0, 0.1)
    
    Hmax.beta[l] ~ dnorm(0, 0.1)
    
  }
  
  
  #################################################################
  ## 9. BETWEEN-CELL VARIANCE
  #################################################################
  
  gamma.tau ~ dgamma(0.01, 0.01)
  
  alpha.tau ~ dgamma(0.01, 0.01)
  
  lambda.ndvi.tau ~ dgamma(0.01, 0.01)
  
  A.tau ~ dgamma(0.01, 0.01)
  
  H0.tau ~ dgamma(0.01, 0.01)
  
  Hmax.tau ~ dgamma(0.01, 0.01)
  
  
  #################################################################
  ## 10. OBSERVATION ERROR
  #################################################################
  
  tau.ndvi ~ dgamma(0.01, 0.01)
  
  tau.height ~ dgamma(0.01, 0.01)
  
  
  #################################################################
  ## 11. STANDARD DEVIATIONS
  #################################################################
  
  sigma.ndvi <- 1 / sqrt(tau.ndvi)
  
  sigma.height <- 1 / sqrt(tau.height)
  
  gamma.sigma <- 1 / sqrt(gamma.tau)
  
  alpha.sigma <- 1 / sqrt(alpha.tau)
  
  lambda.ndvi.sigma <- 1 / sqrt(lambda.ndvi.tau)
  
  A.sigma <- 1 / sqrt(A.tau)
  
  H0.sigma <- 1 / sqrt(H0.tau)
  
  Hmax.sigma <- 1 / sqrt(Hmax.tau)
  
  delta.sigma <- 1 / sqrt(delta.tau)
  
  
  #################################################################
  ## 12. OPTIONAL PREDICTIONS OF HEIGHT THROUGH TIME
  #################################################################
  
  for (i in 1:nGrid) {
    
    for (t in 1:nPredAge) {
      
      recovery.height.pred[i,t] <-
        1 - exp(
          -pred.age[t] / lambda.height[i]
        )
      
      height.pred[i,t] <-
        H0[i] +
        (Hmax[i] - H0[i]) *
        recovery.height.pred[i,t]
      
    }
    
  }
  
}
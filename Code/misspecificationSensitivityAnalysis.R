# Appendix G
# Sensitivity analysis to test robustness to different types of misspecifications
# Figures 9 and 10

library(MASS)
library(dplyr)
library(survival)
library(stringr)
library("mvnTest")
library(gnorm)
library(ggplot2)
library(ggpubr)

######################################################################
# Functions to simulate trial design

# Fit the copula parameters so that our PFS and OS correlation is as
# close as possible to rhoI.
fit.copula <- function(exp.OS.info, ctl.OS.info, exp.PFS.info,
                       ctl.PFS.info, params) {
  grid <- seq(-0.95, 0.95, .01)
  cors <- t(sapply(grid, function(rho) {
    times <- sim.patients(exp.OS.info, ctl.OS.info, exp.PFS.info,
                          ctl.PFS.info, c("rho.exp"=rho, "rho.ctl"=rho),
                          5000, params)
    c("rho.exp"=cor(times$OS[times$isexp == 1], times$PFS[times$isexp == 1]),
      "rho.ctl"=cor(times$OS[times$isexp == 0], times$PFS[times$isexp == 0]))
  }))
  c("rho.exp"=grid[which.min(abs(cors[,"rho.exp"] - params$rhoI))],
    "rho.ctl"=grid[which.min(abs(cors[,"rho.ctl"] - params$rhoI))])
}


# Parameters for our simulation. We start with
# all of them except our study-level prior, which we
# will fit in a short while. We are also missing the
# parameters for type A/B/C, which we will optimize later.
params <- list(lE = 8,
               rhoI = 0.688,
               cp = 50000,
               cw = 20340000,
               alpha = 0.05,
               beta = 0.2,
               delta = 0.5,
               surr.rate = (log(2))/7.2,
               true.rate = (log(2))/21.3,
               gamma = 0.08/12) #check if imputeType actually needed?



# Trial stopping parameters (from section 4, mainSim.R). The results are
# hard-coded here. We use typeC MBC trial parameters assuming no misspecification
params$s0S <- 0.3264536
params$s0T <- 0.2546030
params$rho0 <- 0.7373159
params$n <- 441 
params$v = 2.135741e-02 #no intermediate analysis
params$z = 1.959964 #no intermediate analysis






#Preliminary functions
#What is the expected proportion of patients who have had a surrogate outcome by time t?
E <- function(t, l, lE, n) {
  ret <- lE/(l*n) * ifelse(t <= n/lE, l*t + exp(-l*t) - 1, l*n/lE + exp(-l*t)-exp(-l*(t-n/lE))) 
  ret[is.nan(ret)] <- 0
  ret
}
A <- function(t, muo, lo, lE, n) (E(t, lo, lE, n) + E(t, exp(muo)*lo, lE, n)) / 2

# What is the variance of a type C trial, as a function of its
# number of observations?
vC <- function(nS, nT, s0S, s0T, rho0, rhoI){ 
  D0 <- s0S^2*s0T^2*(1-rho0^2)
  (s0T^2*(1-rhoI^2)+nT*D0/4+(nS-nT)*D0*(1-rhoI^2)/4) / (1-rhoI^2+nT^2*D0/16+nT*(s0T^2+s0S^2-2*rho0*rhoI*s0S*s0T)/4+(nS-nT)*(nT*D0/4+s0S^2*(1-rhoI^2))/4)
}

# What is the expected timing of an analysis that occurs at variance v 
# given that we know muS, muT
tC <- function(muS, muT, v, params) {
  if (params$s0T^2 < v) {
    0  # We immediately hit the number of people needed
  } else {
    tMaxes <- exp(seq(0, log(.Machine$double.xmax), length.out=100))
    tVars <- vC(params$n*A(tMaxes, muS, params$surr.rate, params$lE, params$n), 
                params$n*A(tMaxes, muT, params$true.rate, params$lE, params$n), params$s0S, params$s0T, params$rho0, params$rhoI)
    if (all(tVars >= v)) {
      if (abs(muT) < 100) {
        print(paste("Warning: tC hit max t with muT=", muT))
      }
      return(10^300)
    }
    tMax <- min(tMaxes[tVars < v])
    tryCatch({
      uniroot(function(tt) {
        vC(params$n*A(tt, muS, params$surr.rate, params$lE, params$n), params$n*A(tt, muT, params$true.rate, params$lE, params$n), 
           params$s0S, params$s0T, params$rho0, params$rhoI)-v
      }, c(0, tMax), extendInt = "down", maxiter=10000)$root
    }, error = function(e) { print("C") ; browser() })
  }
}





# returns the indices of the n entries in vals that are closest to target
closest.to <- function(vals, n, target){
  comp <- sapply(vals, function(x){abs(x-target)})
  sort(comp, index.return = TRUE)$ix[1:n]
}







##############################################################################
# helper functions for simulating misspecified trials

#Drawwing effect sizes
draw.effects <- function(num, muS, muT, time, params){
  SI <- cbind(c(4, 4*params$rhoI), c(4*params$rhoI, 4)) 
  
  # Need pJ (proportion of patients w both surrogate & true outcome)
  # and pN (proportion of patients w just surrogate outcome)
  pS <- A(time, muS, params$surr.rate, params$lE, params$n)
  pT <- A(time, muT, params$true.rate, params$lE, params$n)
  qJ <- params$n*(pT)
  qN <- params$n*(pS-pT)
  
  #Draw effect size estimates for trial:
  eJ <- mvrnorm(num, c(muS, muT), SI/(qJ))
  eN <- rnorm(num, muS, sqrt(4/(qN)))
  if(num == 1){
    names(eJ) <- c("eJS", "eJT")
    return(c(eJ, "eN" = eN, "qJ" = qJ, "qN" = qN))
  } else {
    colnames(eJ) <- c("eJS", "eJT")
    return(cbind(eJ, eN, "qJ" = rep(qJ, num), "qN" = rep(qN, num)))
  }
  
  
}


trial.update <- function(ehat, params){
  hat.eJ <- c(ehat["eJS"],ehat["eJT"])
  hat.eN <- c(ehat["eN"], 0)
  S <- rbind(c(1/4, 0), c(0, 0))
  U <- 1/(4*(1-params$rhoI^2))* rbind(c(params$rhoI^2, -1*params$rhoI), c(-1*params$rhoI, 1))
  Sig0 <- with(params, rbind(c(s0S^2, rho0*s0S*s0T), c(rho0*s0S*s0T, s0T^2)))
  final.Sig <- ginv(ginv(Sig0) + (ehat["qJ"]+ ehat["qN"])*S + ehat["qJ"]*U)
  final.mu <- as.vector(final.Sig %*% ((ehat["qJ"])*(S+U) %*% hat.eJ  + ehat["qN"]*S%*%hat.eN))
  
  # Log the final decision, costs, and early stopping.
  final.v <- final.Sig[2,2]
  final.z <- final.mu[2] / sqrt(final.v-(final.v/params$s0T)^2)
  if (!is.finite(final.z)) browser()
  if (abs(final.z) > params$z) {
    reject <- 1
  } else {
    reject <- 0
  }
  
  # Return
  c(c(hat.eJ, hat.eN,
    "final.v"=final.v, "final.z"=final.z, 
      "reject"=reject
  ))
}




# Simulating a typical trial
#   mu: vector of effect size [muS, muT]
#   params: study-level parameters
#   dS: linear shift to surrogate effect size estimate (used when effect size center misspecified)
#   dT: linear shift to effect size estimate (used when effect size center misspecified)
sim.trial.new <- function(mu, params, dS = 0, dT = 0){
  muS <- mu[1]
  muT <- mu[2]
  
  time <- tC(muS, muT, params$v, params)
  
  ehat <- draw.effects(1, muS+dS, muT+dT, time, params)
  
  #outcomes <- t(sapply(1:1000, function(i){
  #  trial.update(ehat[i,], params)
  #}))
  #c(colMeans(outcomes), "muS" = muS, "muT" = muT)
  
  c(trial.update(ehat, params), "muS" = muS, "muT" = muT)
  
}


# Simulating one trial where the prior rho0 is mis-specified
#   mu: vector of effect size [muS, muT]
#   params: study-level parameters
#   rho0new: actual rho0, not misspecified one used for trial design
#   dS: linear shift to surrogate effect size estimate (always 0)
#   dT: linear shift to effect size estimate (always 0)
sim.trial.rho0 <- function(mu, params, rho0new, dS = 0, dT = 0){
  muS <- mu[1]
  muT <- mu[2]
  
  rho0.old <- params$rho0
  params$rho0 <- rho0new
  
  time <- tC(muS, muT, params$v, params)
  
  ehat <- draw.effects(1, muS+dS, muT+dT, time, params)
  
  params$rho0 <- rho0.old
  
  #outcomes <- t(sapply(1:1000, function(i){
  #  trial.update(ehat[i,], params)
  #}))
  #c(colMeans(outcomes), "muS" = muS, "muT" = muT)
  
  c(trial.update(ehat, params), "muS" = muS, "muT" = muT)
  
}


#Simulating a trial where the prior s0S is mis-specified
#   mu: vector of effect size [muS, muT]
#   params: study-level parameters
#   s0Snew: actual s0S, not misspecified one used for trial design
#   dS: linear shift to surrogate effect size estimate (always 0)
#   dT: linear shift to effect size estimate (always 0)
sim.trial.s0S <- function(mu, params, s0Snew, dS = 0, dT = 0){
  muS <- mu[1]
  muT <- mu[2]
  
  s0S.old <- params$s0S
  params$s0S <- s0Snew
  
  time <- tC(muS, muT, params$v, params)
  
  ehat <- draw.effects(1, muS+dS, muT+dT, time, params)
  
  params$s0S <- s0S.old
  
  #outcomes <- t(sapply(1:1000, function(i){
  #  trial.update(ehat[i,], params)
  #}))
  #c(colMeans(outcomes), "muS" = muS, "muT" = muT)
  
  c(trial.update(ehat, params), "muS" = muS, "muT" = muT)
  
}





#Simulate multiple trials with effects drawn from a gaussian prior
set.seed(7)
sim.gaussian <- function(params){
  # Draw effect sizes from Gaussian prior distribution:
    
  #randomly simulate each trial in set A (muT = 0)
  setA <- rnorm(50000, 0, sqrt((1-params$rho0^2)*params$s0S^2))
  A_results <- t(sapply(setA, function(a){
    effect <- c(a,0)
    sim.trial.new(effect, params)
  }))
  
  #randomly simulate each trial in set B (muT = -0.5)
  setB <- rnorm(50000, -1*params$s0S*params$rho0*params$delta/params$s0T, sqrt((1-params$rho0^2)*params$s0S^2))
  B_results <- t(sapply(setB, function(b){
    effect <- c(b, -1*params$delta)
    sim.trial.new(effect, params)
  }))
  
  #randomly simulate each trial in set C
  #setC <- mvrnorm(1000, mu = c(0, 0), Sigma = rbind(c(params$s0S^2, params$rho0*params$s0S*params$s0T), c(params$rho0*params$s0S*params$s0T, params$s0T^2)))
  #C_results <- sapply(1:1000, function(x){sim.trial.new(setC[x,], params)})
  
  typeI <- mean(A_results[,"reject"])
  typeII <- mean(B_results[,"reject"])
  
  c("typeI" = typeI, "typeII" = typeII)
}




# Simulate multiple trials with effects drawn from a gaussian prior 
#   & mis-specified s0S
set.seed(7)
sim.gaussian.s0S <- function(params, s0SNew){
  # Draw effect sizes from Gaussian prior distribution:
  
  #randomly simulate each trial in set A (muT = 0)
  setA <- rnorm(50000, 0, sqrt((1-params$rho0^2)*s0SNew^2))
  A_results <- t(sapply(setA, function(a){
    effect <- c(a,0)
    sim.trial.s0S(effect, params, s0SNew)
  }))
  
  #randomly simulate each trial in set B (muT = -0.5)
  setB <- rnorm(50000, -1*s0SNew*params$rho0*params$delta/params$s0T, sqrt((1-params$rho0^2)*s0SNew^2))
  B_results <- t(sapply(setB, function(b){
    effect <- c(b, -1*params$delta)
    sim.trial.s0S(effect, params, s0SNew)
  }))
  
  #randomly simulate each trial in set C
  #setC <- mvrnorm(1000, mu = c(0, 0), Sigma = rbind(c(params$s0S^2, params$rho0*params$s0S*params$s0T), c(params$rho0*params$s0S*params$s0T, params$s0T^2)))
  #C_results <- sapply(1:1000, function(x){sim.trial.new(setC[x,], params)})
  
  typeI <- mean(A_results[,"reject"])
  typeII <- mean(B_results[,"reject"])
  
  c("typeI" = typeI, "typeII" = typeII)
}


# Simulate multiple trials with effects drawn from a gaussian prior 
#   & misspecified rho0
set.seed(7)
sim.gaussian.rho0 <- function(params, rho0New){
  # Draw effect sizes from Gaussian prior distribution:
  
  #randomly simulate each trial in set A (muT = 0)
  setA <- rnorm(50000, 0, sqrt((1-rho0New^2)*params$s0S^2))
  A_results <- t(sapply(setA, function(a){
    effect <- c(a,0)
    sim.trial.rho0(effect, params, rho0New)
  }))
  
  #randomly simulate each trial in set B (muT = -0.5)
  setB <- rnorm(50000, -1*params$s0S*rho0New*params$delta/params$s0T, sqrt((1-rho0New^2)*params$s0S^2))
  B_results <- t(sapply(setB, function(b){
    effect <- c(b, -1*params$delta)
    sim.trial.rho0(effect, params, rho0New)
  }))
  
  #randomly simulate each trial in set C
  #setC <- mvrnorm(1000, mu = c(0, 0), Sigma = rbind(c(params$s0S^2, params$rho0*params$s0S*params$s0T), c(params$rho0*params$s0S*params$s0T, params$s0T^2)))
  #C_results <- sapply(1:1000, function(x){sim.trial.new(setC[x,], params)})
  
  typeI <- mean(A_results[,"reject"])
  typeII <- mean(B_results[,"reject"])
  
  c("typeI" = typeI, "typeII" = typeII)
}



# Simulate multiple trials with effects drawn from a gaussian prior 
#   and with mean surrogate (true) outcome effect sizes shifted over by dS (dT)
set.seed(7)
sim.gaussian.effectShift <- function(params, dS, dT){
  # Draw effect sizes from Gaussian prior distribution:
  
  #randomly simulate each trial in set A (muT = 0)
  setA <- rnorm(50000, 0, sqrt((1-params$rho0^2)*params$s0S^2))
  A_results <- t(sapply(setA, function(a){
    effect <- c(a,0)
    sim.trial.new(effect, params, dS, dT)
  }))
  
  #randomly simulate each trial in set B (muT = -0.5)
  setB <- rnorm(50000, -1*params$s0S*params$rho0*params$delta/params$s0T, sqrt((1-params$rho0^2)*params$s0S^2))
  B_results <- t(sapply(setB, function(b){
    effect <- c(b, -1*params$delta)
    sim.trial.new(effect, params, dS, dT)
  }))
  
  #randomly simulate each trial in set C
  #setC <- mvrnorm(1000, mu = c(0, 0), Sigma = rbind(c(params$s0S^2, params$rho0*params$s0S*params$s0T), c(params$rho0*params$s0S*params$s0T, params$s0T^2)))
  #C_results <- sapply(1:1000, function(x){sim.trial.new(setC[x,], params)})
  
  typeI <- mean(A_results[,"reject"])
  typeII <- mean(B_results[,"reject"])
  
  c("typeI" = typeI, "typeII" = typeII, "dS" = dS, "dT" = dT)
}




# Simulate multiple trials with a effect sizes drawn from a generalized 
#   normal distribution
set.seed(7)
# muS chosen from general normal distribution
sim.general <- function(params, shape, mean.change = 0){
  # Draw muS from generalized normal prior distribution:
  
  #randomly simulate each trial in set A (muT = 0)
  alphaS <- sqrt(2)*sqrt((1-params$rho0^2)*params$s0S^2)
  setA <- rgnorm(50000, mu= mean.change, alpha = alphaS, beta = shape)
  A_results <- t(sapply(setA, function(a){
    effect <- c(a,0)
    sim.trial.new(effect, params)
  }))
  
  #randomly simulate each trial in set B (muT = -0.5)
  alphaS <- sqrt(2)*sqrt((1-params$rho0^2)*params$s0S^2)
  setB <- rgnorm(50000, mu = mean.change+-1*params$s0S*params$rho0*params$delta/params$s0T, 
                alpha = alphaS, beta = shape)
  B_results <- t(sapply(setB, function(b){
    effect <- c(b, -1*params$delta)
    sim.trial.new(effect, params)
  }))
  
  #randomly simulate each trial in set C
  #setC <- mvrnorm(1000, mu = c(0, 0), Sigma = rbind(c(params$s0S^2, params$rho0*params$s0S*params$s0T), c(params$rho0*params$s0S*params$s0T, params$s0T^2)))
  #C_results <- sapply(1:1000, function(x){sim.trial.new(setC[x,], params)})
  
  typeI <- mean(A_results[,"reject"])
  typeII <- mean(B_results[,"reject"])
  
  c("typeI" = typeI, "typeII" = typeII, "shape" = shape)
}








###############################################################################
# Actually run simulations of type C trials with various misspecifications 
#   create subfigures for Figure 9, create 1 subfigure for figure 10
choices <- seq(0.9, 1.1, by = 0.01)
choices_bias <- seq(-0.5, 0.5, by = 0.05)

# Actually simulate with effect sizes pulled from generalized normal distribution
#   with changing shape parameter. (remember that regular normal is generalized
#   normal with shape parameter 2 & scale parameter alpha^2 = 2*sigma^2)
gen.shapeChange <- t(sapply(choices, function(x){
  shape <- 2*x
  
  results <- c(sim.general(params, shape), "percent" = x)
}))

shape_plot <- ggplot(as.data.frame(gen.shapeChange), aes(100*(percent-1), typeI)) + geom_point() + geom_smooth(aes(colour="Type I"))+ 
  geom_point(aes(100*(percent-1), 1-typeII)) + geom_smooth(aes(100*(percent-1), 1-typeII, colour = "Type II")) +
  geom_vline(xintercept = 0) +
  labs(x = "Change to shape parameter (%)", y = "Error rate")


#Actually simulate with changing mean:
norm.meanChange <- t(sapply(choices_bias, function(x){
  mean.change <- x*params$s0S
  results <- c(sim.general(params, 2, mean.change), "mean.change" = mean.change, "percent" = x)
}))

mean_plot <- ggplot(as.data.frame(norm.meanChange), aes(percent, typeI)) + geom_point() + geom_smooth(aes(colour="Type I"))+ 
  geom_point(aes(percent, 1-typeII)) + geom_smooth(aes(percent, 1-typeII, colour = "Type II")) +
  geom_vline(xintercept = 0) +
  labs(x = expression(paste("Bias added to ", mu[S], " (units of ",  sigma['0S']*'', ")")), y = "Error rate")



#Actually simulate with changing variance:
s0S.orig <- params$s0S
norm.varChange <- t(sapply(choices, function(x){
  News0S <- s0S.orig*(x)
  print(News0S)
  results <- c(sim.gaussian.s0S(params, News0S), "s0S" = News0S, "percent" = x)
}))

params$s0S <- s0S.orig

var_plot <- ggplot(as.data.frame(norm.varChange), aes(100*(percent-1), typeI)) + geom_point() + geom_smooth(aes(colour = "Type I")) +
  geom_point(aes(100*(percent-1), 1-typeII)) + geom_smooth(aes(100*(percent-1), 1-typeII, colour = "Type II")) +
  geom_vline(xintercept = 0) +
  labs(x = expression(paste("Change to ", sigma['0S']*'', " (%)")), y = "Error rate")



# Actually simulate with changing underlying rho0
choices2 <- seq(0.95, 1.05, by = 0.005)
rho0.orig <- params$rho0
norm.rho0Change <- t(sapply(choices, function(x){
  rho0new <- rho0.orig*x
  results <- c(sim.gaussian.rho0(params, rho0new), "rho0" = rho0new, "percent" = x)
}))
params$rho0 <- rho0.orig

rho0_plot <- ggplot(as.data.frame(norm.rho0Change), aes(100*(percent-1), typeI)) + geom_point() + geom_smooth(aes(colour = "Type I")) +
  geom_point(aes(100*(percent-1), 1-typeII)) + geom_smooth(aes(100*(percent-1), 1-typeII, colour = "Type II")) +
  geom_vline(xintercept = 0) +
  labs(x = expression(paste("Change to ", rho[0], " (%)")), y = "Error rate")



# Actually simulate with varying underlying rhoI
rhoI.orig <- params$rhoI
norm.rhoIChange <- t(sapply(choices2, function(x){
  params$rhoI <- rhoI.orig*x
  results <- c(sim.gaussian(params), "rhoI" = params$rhoI, "percent" = x)
}))
params$rhoI <- rhoI.orig

rhoI_plot <- ggplot(as.data.frame(norm.rhoIChange), aes(100*(1-percent), typeI)) + geom_point() + geom_smooth(aes(colour = "Type I")) +
  geom_point(aes(100*(1-percent), 1-typeII)) + geom_smooth(aes(100*(1-percent), 1-typeII, colour = "Type II")) +
  geom_vline(xintercept = 0) +
  labs(x = expression(paste("Change to ", rho[I], " (%)")), y = "Error rate")




#Actually simulate with changing effect size center.

#   Shift surrogate effect size center:
norm.effectShiftS <- t(sapply(choices_bias, function(x){
  delta <- x*params$s0S
  results <- c(sim.gaussian.effectShift(params, delta, 0), "percentage" = x)
}))

effBiasS_plot <- ggplot(as.data.frame(norm.effectShiftS), aes(percentage, typeI)) + geom_point() + geom_smooth(aes(colour = "Type I")) + 
  geom_point(aes(percentage, 1-typeII)) + geom_smooth(aes(percentage, 1-typeII, colour = "Type II")) + 
  geom_vline(xintercept =0) +
  labs(x = expression(paste("Bias added to ", e[S], " (units of ",  sigma['0S']*'', ")")), y = "Error rate")

#   Shift True outcome effect size center
norm.effectShiftT <- t(sapply(choices_bias, function(x){
  delta <- x*params$s0T
  results <- c(sim.gaussian.effectShift(params, 0, delta), "percentage" = x)
}))

effBiasT_plot <- ggplot(as.data.frame(norm.effectShiftT), aes(percentage, typeI)) + geom_point() + geom_smooth(aes(colour = "Type I")) + 
  geom_point(aes(percentage, 1-typeII)) + geom_smooth(aes(percentage, 1-typeII, colour = "Type II")) + 
  geom_vline(xintercept =0) +
  labs(x = expression(paste("Bias added to ", e[T], " (units of ",  sigma['0T']*'', ")")), y = "Error rate") + ylim(0, 0.45)








# Figure 9
pdf("../Plots/misspecificationSensitivity.pdf", width = 10, height=10)
ggarrange(rhoI_plot, rho0_plot, var_plot, shape_plot, mean_plot, effBiasS_plot, 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3, legend = "none")
dev.off()











#################################################################################
# Optimal parameters for a Type A trial with no misspecification

paramsA <- params
paramsA$n <- 424
paramsA$v <- 2.135741e-02
paramsA$z <- 1.959964



#################################################################################
# Helper functions for simulating mis-specified Type A trial

tA <- function(muT, v, n, p) {
  if (p$s0T^2 < v) {
    0  # We immediately hit the number of people needed
  } else {
    tMaxes <- exp(seq(0, log(.Machine$double.xmax), length.out=100))
    tVars <- vA(n*A(tMaxes, muT, p$true.rate, p$lE, n), p$s0T)
    if (all(tVars >= v)) {
      browser()  # Should not happen -- none of our designs hit variance
    }
    tMax <- min(tMaxes[tVars < v])
    tryCatch({
      uniroot(function(tt) {
        vA(n*A(tt, muT, p$true.rate, p$lE, n), p$s0T)-v
      }, c(0, tMax), extendInt = "down", maxiter=10000)$root
    }, error = function(e) { print("A") ; browser() })
  }
}
vA <- function(nT, s0T) s0T^2/(1+nT*s0T^2/4)

trial.update.A <- function(ehat, params){
  hat.eJ <- c(ehat["eJS"],ehat["eJT"])
  hat.eN <- c(ehat["eN"], 0)
  T <- rbind(c(0, 0), c(0, 1/4))
  Sig0 <- with(params, rbind(c(s0S^2, rho0*s0S*s0T), c(rho0*s0S*s0T, s0T^2)))
  final.Sig <- ginv(ginv(Sig0) + (ehat["qJ"])*T)
  final.mu <- as.vector(final.Sig %*% ((ehat["qJ"])*(T) %*% hat.eJ))
  
  # Log the final decision, costs, and early stopping.
  final.v <- final.Sig[2,2]
  final.z <- final.mu[2] / sqrt(final.v-(final.v/params$s0T)^2)
  if (!is.finite(final.z)) browser()
  if (abs(final.z) > params$z) {
    reject <- 1
  } else {
    reject <- 0
  }
  
  # Return
  c(c(hat.eJ, hat.eN,
      "final.v"=final.v, "final.z"=final.z, 
      "reject"=reject
  ))
}





#################################################################################
# Functions to simulate a misspecified type A trial

sim.trial.new.A <- function(mu, params, dS = 0, dT = 0){
  muS <- mu[1]
  muT <- mu[2]
  
  time <- tA(muT, params$v, params$n, params)
  
  ehat <- draw.effects(1, muS+dS, muT+dT, time, params)
  
  #outcomes <- t(sapply(1:1000, function(i){
  #  trial.update(ehat[i,], params)
  #}))
  #c(colMeans(outcomes), "muS" = muS, "muT" = muT)
  
  c(trial.update.A(ehat, params), "muS" = muS, "muT" = muT)
  
}


set.seed(7)
sim.gaussian.effectShift.A <- function(params, dS, dT){
  # Draw effect sizes from Gaussian prior distribution:
  
  #randomly simulate each trial in set A (muT = 0)
  setA <- rnorm(50000, 0, sqrt((1-params$rho0^2)*params$s0S^2))
  A_results <- t(sapply(setA, function(a){
    effect <- c(a,0)
    sim.trial.new.A(effect, params, dS, dT)
  }))
  
  #randomly simulate each trial in set B (muT = -0.5)
  setB <- rnorm(50000, -1*params$s0S*params$rho0*params$delta/params$s0T, sqrt((1-params$rho0^2)*params$s0S^2))
  B_results <- t(sapply(setB, function(b){
    effect <- c(b, -1*params$delta)
    sim.trial.new.A(effect, params, dS, dT)
  }))
  
  #randomly simulate each trial in set C
  #setC <- mvrnorm(1000, mu = c(0, 0), Sigma = rbind(c(params$s0S^2, params$rho0*params$s0S*params$s0T), c(params$rho0*params$s0S*params$s0T, params$s0T^2)))
  #C_results <- sapply(1:1000, function(x){sim.trial.new(setC[x,], params)})
  
  typeI <- mean(A_results[,"reject"])
  typeII <- mean(B_results[,"reject"])
  
  c("typeI" = typeI, "typeII" = typeII, "dS" = dS, "dT" = dT)
}



###############################################################################
# Actually run simulations of type C trials with various misspecifications 
#   create subfigures for Figure 9, create 1 subfigure for figure 10

# Actually simulate with changing effect size center.
#   Shift true outcome effect size center
norm.effectShiftT.A <- t(sapply(choices_bias, function(x){
  delta <- x*paramsA$s0T
  results <- c(sim.gaussian.effectShift.A(paramsA, 0, delta), "percentage" = x)
}))

effBiasT_plot.A <- ggplot(as.data.frame(norm.effectShiftT.A), aes(percentage, typeI)) + geom_point() + geom_smooth(aes(colour = "Type I")) + 
  geom_point(aes(percentage, 1-typeII)) + geom_smooth(aes(percentage, 1-typeII, colour = "Type II")) + 
  geom_vline(xintercept =0) +
  labs(x = expression(paste("Bias added to ", e[T], " (units of ",  sigma['0T']*'', ")")), y = "Error rate") + ylim(0, 0.45)



# Figure 10
          
pdf("../Plots/misspecificationSensitivity_effBias.pdf", width = 10, height=4)
ggarrange(effBiasT_plot.A, effBiasT_plot, 
          labels = c("Type A Trial", "Type C Trial"),
          ncol = 2, nrow = 1, legend = "none")
dev.off()



#Save data

write.csv(norm.rhoIChange, "../Simulations/rhoIMisspec.csv")
write.csv(norm.rho0Change, "../Simulations/rho0Misspec.csv")
write.csv(norm.meanChange, "../Simulations/meanMisspec.csv")
write.csv(norm.varChange, "../Simulations/varMisspec.csv")
write.csv(gen.shapeChange, "../Simulations/distMisspec.csv")
write.csv(norm.effectShiftS, "../Simulations/effBiasS.csv")


write.csv(norm.effectShiftT, "../Simulations/effBiasT.csv")
write.csv(norm.effectShiftT.A, "../Simulations/effBiasT_A.csv")



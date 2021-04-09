# Code for Section 4, main simulation using actual trial data
# Simulation results
# Table 4.3
# Figures 5 & 6

library("XML")
library(MASS)
library(dplyr)
library(survival)
library(stringr)
library("mvnTest")
library(ggplot2)
setwd("../Data/MBC Data")

#Helper functions

#Trim filenames
trim.trailing <- function (x) sub("\\s+$", "", x)
##############################################
# Load the xml of an extracted K-M curve, clean up inconsistencies, and then return
# both an interpolation function as well as the maximum follow-up time.
extract.km <- function(fname) {
  # Grab the points off the K-M curve (obtained with PlotDigitizer java app)
  parsed <- xmlToList(xmlParse(fname))
  points <- do.call(rbind, parsed[names(parsed) == "point"])[,c("dx", "dy")]
  mode(points) <- "numeric"
  
  # Convert based on the units for the x- and y-axis (see checkKM.R where we
  # dug into all of this; we assume all arms being processed have passed the checks
  # in that script file, so there are no data validity checks in this simulation script)
  an <- str_trim(tolower(parsed$axesnames))
  if (an[1] %in% c("days", "time (in days)")) {
    points[,1] <- points[,1] / 30
  } else if (an[1] %in% c("years", "time from random assignment (years)", "years from randomization")) {
    points[,1] <- points[,1] * 12
  } else if (an[1] %in% c("weeks", "week")) {
    points[,1] <- points[,1] * 7 / 30
  }
  if (an[2] %in% c("percent", "%", "pfs (%)", "overall survival (%)", "percent of progression-free") || (max(points[,2]) > 90 && max(points[,2]) < 110)) {
    points[,2] <- points[,2] / 100
  }
  
  # Clean up data inconsistencies
  points[,1] <- pmax(points[,1], 0)
  points[,2] <- pmax(pmin(points[,2], 1), 0)
  points <- points[points[,1] > 0,]
  points <- rbind(c(0, 1), points)
  
  # If y increases in consecutive points A and B, then increase A's y to match B
  while (any(diff(points[,2]) > 0)) {
    inc.diff <- diff(points[,2]) > 0
    points[c(inc.diff, FALSE),2] <- points[c(FALSE, inc.diff),2]
  }
  
  # If x decreases in consecutive points A and B, then decrease A's x to match B
  while (any(diff(points[,1]) < 0)) {
    inc.diff <- diff(points[,1]) < 0
    points[c(inc.diff, FALSE),1] <- points[c(FALSE, inc.diff),1]
  }
  
  
  
  list(interp = approxfun(points[,2], points[,1]), points = points,
       max.followup = max(points[,1]),
       fname = fname)
}

##############################################
#Load MBC data and clean up a bit (remove non-randomized studies or randomized studies with only a single
# arm present; it sounds impossible to be randomized but only have a single arm, but sometimes our
# exclusion factors match some but not all arms, leaving one left from a RCT)
breast.raw <- read.csv("MBC_data_final.csv", stringsAsFactors = FALSE)
# can use either TTP or PFS as surrogate, since similar measurements
breast.raw$surr_KM <- ifelse(breast.raw$PFS_KM != "", breast.raw$PFS_KM, breast.raw$TTP_KM)
breast.raw <- subset(breast.raw, Randomized == "Yes")
breast.phase <- read.csv("MBC_phase.csv", stringsAsFactors = FALSE)
breast.raw <- merge(breast.raw, breast.phase, by="Unique_ID")
# rcode identifies study and randomization group (some studies randomize A vs. B and C vs. D, so there are
# 2 groups separately being randomized)
breast.raw$rcode <- paste(trim.trailing(as.character(sapply(strsplit(as.character(breast.raw$Unique_ID), "_"), "[", 1))), breast.raw$Rand_Group)

#Make sure that we correctly define which arms are control & which experiemnt
flip.codes <- function(x, dat) {
  exptype <- dat$Arm_Type[match(x$IDexp, dat$Unique_ID)]
  ctltype <- dat$Arm_Type[match(x$IDctl, dat$Unique_ID)]
  to.flip <- exptype == "control"
  flip.exp <- x$IDexp[to.flip]
  x$IDexp[to.flip] <- x$IDctl[to.flip]
  x$IDctl[to.flip] <- flip.exp
  x$effect[to.flip] <- -x$effect[to.flip]
  x$armcode <- paste(x$IDexp, "|", x$IDctl)
  x
}

# Summarize the data for each rcode and figure out which rcodes meet all our
# criteria for inclusion (have an OS and surrogate KM curve, have exactly 2 arms with exactly
# one control and one experiment, and have a reported HR for OS and surrogate)
summary.rcode <- do.call(rbind, lapply(split(breast.raw, breast.raw$rcode), function(x) {
  data.frame(rcode=x$rcode[1],
             arms=nrow(x),
             numOSKM=sum(x$OS_KM != ""),
             numsurrKM=sum(x$surr_KM != ""),
             phase=x$Phase[1],
             numControl=sum(x$Arm_Type == "control"),
             numExp=sum(x$Arm_Type == "experiment"),
             numUnlabeled=sum(x$Arm_Type != "control" & x$Arm_Type != "experiment"),
             stringsAsFactors=FALSE)
}))
summary.rcode <- subset(summary.rcode, arms == 2 & numOSKM == 2 & numsurrKM == 2 &
                          numControl == 1 & numExp == 1)
keep.rcode <- summary.rcode$rcode



####################################################################################
# Actually Load the data

# Load the data. This object holds interpolation functions for each KM curve from 
#   the trials we want to simulate; will use this as the empirical distribution.
big.res <- lapply(keep.rcode, function(x) {
  ctl.info <- subset(breast.raw, rcode == x & Arm_Type == "control")
  exp.info <- subset(breast.raw, rcode == x & Arm_Type == "experiment")
  
  # Grab the KM curves for the control and experiment arms
  error <- FALSE
  tryCatch({
    ctl.OS.KM <- extract.km(paste0("KM/", ctl.info$OS_KM))
    ctl.surr.KM <- extract.km(paste0("KM/", ctl.info$surr_KM))
    exp.OS.KM <- extract.km(paste0("KM/", exp.info$OS_KM))
    exp.surr.KM <- extract.km(paste0("KM/", exp.info$surr_KM))
  }, error = function(e) {
    print(e)
    print(paste("Error extracting KM curve in", ctl.info$rcode))
    error <<- TRUE
  })
  if (error) return(NULL)
  list(ctl.OS.KM=ctl.OS.KM, ctl.surr.KM=ctl.surr.KM,
       exp.OS.KM=exp.OS.KM, exp.surr.KM=exp.surr.KM)
})
names(big.res) <- keep.rcode
big.res <- big.res[!sapply(big.res, is.null)]



####################################################################
# Simulation Helper functions

# OK -- we loaded the KM curve and now we want to simulate from
# it. The key challenge is how to deal with the fact that the KM
# curve often ends before it reaches a y-axis value of 0. Let's say
# it ends at a 0.2 value. The question is how do we simulate this
# remaining 20% of the population? 
# 1) for the control and treatment arm, compute the maximum follow-up,
#    and then take the smaller of those two values. Make this be a
#    hard cutoff, where anything in either arm beyond that value is
#    chopped off to that value.
# Now we will have a true event time for every
# patient. I would argue this is an important feature of any
# simulation model -- we need to know the ground truth!
# 
# We also apply our Gaussian copula to link PFS and OS, using rejection
# sampling to reject anything with PFS > OS. This will cause us to not
# exactly match the marginal distributions, but for large copula rho
# it should be pretty close.
# 
# Also randomly draw an enrollment time, using params$lE as the rate
# per time unit (month).
# 
# "expand" is a hidden parameter that handles how much we should expand
# N beyond its requested size. It will be scaled up when needed

# Simulate patients from KM curves
#   exp.OS.info: empirical distribution for OS in experiment arm
#   ctl.OS.info: empirical distribution for OS in control arm
#   exp.PFS.info: empirical distribution for PFS (or TTP) in experiment arm
#   ctl.PFS.info: empirical distribution for PFS (or TTP) in control arm
#   copula: the rho values to ensure the proper correlation between OS & PFS effects
#   N: The target number of patients we want to simulate
#   params: underlying parameters
sim.patients <- function(exp.OS.info, ctl.OS.info, exp.PFS.info,
                         ctl.PFS.info, copula, N, params,
                         expand = 5) {
  # Step 1 -- simulation from the KM curves. We will get "NA" as the
  # time for patients who are simulated to be beyond the maximum
  # follow-up time on the KM curve. Use the Gaussian copula to draw
  # the survival times, and rejection sample to limit to PFS/OS pairs
  # where PFS <= OS.
  cop.exp <- mvrnorm(expand*N, mu = c(0, 0), Sigma = rbind(c(1, copula["rho.exp"]), c(copula["rho.exp"], 1)))
  times.exp <- data.frame(OS = exp.OS.info$interp(pnorm(cop.exp[,2])),
                          PFS = exp.PFS.info$interp(pnorm(cop.exp[,1])),
                          isexp = 1, enroll = runif(expand*N, 0, N/params$lE))
  cop.ctl <- mvrnorm(expand*N, mu = c(0, 0), Sigma = rbind(c(1, copula["rho.ctl"]), c(copula["rho.ctl"], 1)))
  times.ctl <- data.frame(OS = ctl.OS.info$interp(pnorm(cop.ctl[,2])),
                          PFS = ctl.PFS.info$interp(pnorm(cop.ctl[,1])),
                          isexp = 0, enroll = runif(expand*N, 0, N/params$lE))

  # Step 2 -- apply our selected imputation function
  
  cutoff.OS <- min(max(exp.OS.info$points[,"dx"]),
                   max(ctl.OS.info$points[,"dx"]))
  cutoff.PFS <- min(max(exp.PFS.info$points[,"dx"]),
                    max(ctl.PFS.info$points[,"dx"]))
  times <- rbind(times.exp, times.ctl)
  times$OS[is.na(times$OS)] <- cutoff.OS
  times$OS <- pmin(times$OS, cutoff.OS)
  times$PFS[is.na(times$PFS)] <- cutoff.PFS
  times$PFS <- pmin(times$PFS, cutoff.PFS)
  
  
  # Step 3: filter down to those respecting PFS < OS, limit to the
  # requested number, and return
  times <- times[times$PFS <= times$OS,]
  n.exp <- rbinom(1, N, 0.5)  # Number in experimental arm
  if (sum(times$isexp == 1) < n.exp || sum(times$isexp == 0) < N-n.exp) {
    if (expand > 1000) {
      browser() # Huge number with PFS > OS; let's investigate
    }
    sim.patients(exp.OS.info, ctl.OS.info, exp.PFS.info, ctl.PFS.info,
                 copula, N, params, expand*10)
  } else {
    rbind(head(times, n.exp), tail(times, N-n.exp))
  }
}

# Fit the copula parameters so that our PFS and OS correlation is as
# close as possible to rhoI.
#   exp.OS.info: empirical distribution for OS in experiment arm
#   ctl.OS.info: empirical distribution for OS in control arm
#   exp.PFS.info: empirical distribution for PFS (or TTP) in experiment arm
#   ctl.PFS.info: empirical distribution for PFS (or TTP) in control arm
#   params: underlying parameters
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

# Take a "times" object (the OS, PFS, group, and enrollment times) and
# compute the effect sizes through a given time. We will have only
# observed a PFS event if the enrollment time plus the PFS is before
# the analysis time (and similarly for OS).
# returns the hazard ratios, cox proportional hazards, and log-rank
effect.size <- function(times, analysis.time) {
  # Use the analysis time to remove people who are not yet enrolled
  # and to censor those who are enrolled
  times <- times[times$enroll < analysis.time,]
  times$OS.event <- as.numeric(times$enroll + times$OS <= analysis.time)
  times$PFS.event <- as.numeric(times$enroll + times$PFS <= analysis.time)
  times$OS <- times$OS - pmax(0, times$enroll + times$OS - analysis.time)
  times$PFS <- times$PFS - pmax(0, times$enroll + times$PFS - analysis.time)
  
  # Compute the hazard ratios
  S.OS <- Surv(time=times$OS, event=times$OS.event)
  res.OS <- coxph(S.OS~times$isexp)
  S.PFS <- Surv(time=times$PFS, event=times$PFS.event)
  res.PFS <- coxph(S.PFS~times$isexp)
  if (sum(times$OS.event) > 0) {
    S.PFSj <- Surv(time=times$PFS[times$OS.event == 1], event=times$PFS.event[times$OS.event == 1])
    res.PFSj <- coxph(S.PFSj~times$isexp[times$OS.event == 1])
  } else {
    res.PFSj <- list("coefficients"=NA, var=matrix(NA))
  }
  if (sum(1-times$OS.event) > 0) {
    S.PFSs <- Surv(time=times$PFS[times$OS.event == 0], event=times$PFS.event[times$OS.event == 0])
    res.PFSs <- coxph(S.PFSs~times$isexp[times$OS.event == 0])
  } else {
    res.PFSs <- list("coefficients"=NA, var=matrix(NA))
  }
  
  # Logrank test for PFS
  # S.PFS2 <- survdiff(S.PFS~times$isexp)
  # (S.PFS2$obs[2]-S.PFS2$exp[2])/S.PFS2$var[2,2]
  times <- times[order(times$PFS),]
  times$n.exp <- sum(times$isexp) - head(cumsum(c(0, times$isexp)), -1)
  times$n.ctl <- sum(1-times$isexp) - head(cumsum(c(0, 1-times$isexp)), -1)
  O <- with(times, isexp[PFS.event == 1])
  E <- with(times, (n.exp/(n.exp+n.ctl))[PFS.event == 1])
  LR.mu.PFS <- (sum(O)-sum(E)) / sum(E*(1-E))  # HR estimate; matches (S.PFS2$obs[2]-S.PFS2$exp[2])/S.PFS2$var[2,2]
  LR.V.PFS <- 1 / sum(E*(1-E))  # var estimate; matches 1 / S.PFS2$var[2,2]
  if (sum(times$OS.event) > 0) {
    Oj <- with(times, isexp[PFS.event == 1 & OS.event == 1])
    Ej <- with(times, (n.exp/(n.exp+n.ctl))[PFS.event == 1 & OS.event == 1])
    LR.mu.PFSj <- (sum(Oj)-sum(Ej)) / sum(Ej*(1-Ej))
    LR.V.PFSj <- 1 / sum(Ej*(1-Ej))
  } else {
    LR.mu.PFSj <- NA
    LR.V.PFSj <- NA
  }
  if (sum(1-times$OS.event) > 0) {
    Os <- with(times, isexp[PFS.event == 1 & OS.event == 0])
    Es <- with(times, (n.exp/(n.exp+n.ctl))[PFS.event == 1 & OS.event == 0])
    LR.mu.PFSs <- (sum(Os)-sum(Es)) / sum(Es*(1-Es))
    LR.V.PFSs <- 1 / sum(Es*(1-Es))
  } else {
    LR.mu.PFSs <- NA
    LR.V.PFSs <- NA
  }
  
  # Cox Proportion Hazards model for PFS (commented out part reproduces the coxph results above)
  # toopt <- function(b) with(times, sum((n.exp*exp(b)/(n.ctl+n.exp*exp(b)))[PFS.event == 1]))
  # coefficient estimate; matches unname(res.PFS$coefficients)
  # mu.PFS <- uniroot(function(b) toopt(b)-sum(times$PFS.event == 1 & times$isexp == 1), c(-10, 10))$root
  # pexp <- with(times, n.exp*exp(mu.PFS)/(n.ctl+n.exp*exp(mu.PFS)))
  # variance estimate; matches res.PFS$var[1,1]
  # V.PFS <- 1/with(times, sum((pexp*(1-pexp))[PFS.event == 1]))
  if (sum(times$OS.event) > 0) {
    tooptj <- function(b) with(times, sum((n.exp*exp(b)/(n.ctl+n.exp*exp(b)))[PFS.event == 1 & OS.event == 1]))
    C2.mu.PFSj <- uniroot(function(b) tooptj(b)-sum(times$PFS.event == 1 & times$OS.event == 1 & times$isexp == 1), c(-10, 10))$root
    pexpj <- with(times, n.exp*exp(C2.mu.PFSj)/(n.ctl+n.exp*exp(C2.mu.PFSj)))
    C2.V.PFSj <- 1/with(times, sum((pexpj*(1-pexpj))[PFS.event == 1 & OS.event == 1]))
  } else {
    C2.mu.PFSj <- NA
    C2.V.PFSj <- NA
  }
  if (sum(1-times$OS.event) > 0) {
    toopts <- function(b) with(times, sum((n.exp*exp(b)/(n.ctl+n.exp*exp(b)))[PFS.event == 1 & OS.event == 0]))
    C2.mu.PFSs <- uniroot(function(b) toopts(b)-sum(times$PFS.event == 1 & times$OS.event == 0 & times$isexp == 1), c(-10, 10))$root
    pexps <- with(times, n.exp*exp(C2.mu.PFSs)/(n.ctl+n.exp*exp(C2.mu.PFSs)))
    C2.V.PFSs <- 1/with(times, sum((pexps*(1-pexps))[PFS.event == 1 & OS.event == 0]))
  } else {
    C2.mu.PFSs <- NA
    C2.V.PFSs <- NA
  }

  # Return effect sizes and variances
  c("mu.OS"=unname(res.OS$coefficients),
    "V.OS"=res.OS$var[1,1],
    "mu.PFS"=unname(res.PFS$coefficients),
    "V.PFS"=res.PFS$var[1,1],
    "mu.PFSj"=unname(res.PFSj$coefficients),
    "V.PFSj"=res.PFSj$var[1,1],
    "mu.PFSs"=unname(res.PFSs$coefficients),
    "V.PFSs"=res.PFSs$var[1,1],
    "n.OS"=sum(times$OS.event == 1),
    "n.PFS"=sum(times$PFS.event == 1),
    "LR.mu.PFS"=LR.mu.PFS,
    "LR.V.PFS"=LR.V.PFS,
    "LR.mu.PFSj"=LR.mu.PFSj,
    "LR.V.PFSj"=LR.V.PFSj,
    "LR.mu.PFSs"=LR.mu.PFSs,
    "LR.V.PFSs"=LR.V.PFSs,
    "C2.mu.PFSj"=C2.mu.PFSj,
    "C2.V.PFSj"=C2.V.PFSj,
    "C2.mu.PFSs"=C2.mu.PFSs,
    "C2.V.PFSs"=C2.V.PFSs)
}



#######################################################################
# Simulation preliminaries: parameters & arm-specific copulas

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
               gamma = 0.08/12)

# The first step is to learn our copula parameters for the control and
# treatment arms of each study. This relies on parameters rhoI (copula
# gets fit to this correlation), as well as imputeType (used to
# simulate the populations).
set.seed(144)
for (rcode in names(big.res)) {
  print(rcode)
  x <- big.res[[rcode]]
  big.res[[rcode]]$copula <-
    fit.copula(x$exp.OS.KM, x$ctl.OS.KM, x$exp.surr.KM,
               x$ctl.surr.KM, params)
}



########################################################################
#Obtain "ground truth" HR estimates for surrogate and OS for all studies

# Do the master run (100k person trial; assumed ground truth).
# Relies on imputeType only from the parameters.
set.seed(144)
for (rcode in names(big.res)) {
  print(rcode)
  x <- big.res[[rcode]]
  times <- sim.patients(x$exp.OS.KM, x$ctl.OS.KM, x$exp.surr.KM,
                        x$ctl.surr.KM, x$copula, 100000, params)
  big.res[[rcode]]$bigES <- effect.size(times, 1e6)
}

# Do a master run, where we run until a 40-month
# followup. We will scale up the enrollment rate massively
# so that we enroll in the time period [0, 40].
set.seed(144)
old.lE <- params$lE
params$lE <- 100000/40
for (rcode in names(big.res)) {
  print(rcode)
  x <- big.res[[rcode]]
  times <- sim.patients(x$exp.OS.KM, x$ctl.OS.KM, x$exp.surr.KM,
                        x$ctl.surr.KM, x$copula, 100000, params)
  big.res[[rcode]]$bigES40 <- effect.size(times, 40)
}
params$lE <- old.lE


all.sim <- as.data.frame(t(sapply(big.res, "[[", "bigES")))
all.sim40 <- as.data.frame(t(sapply(big.res, "[[", "bigES40")))


# Table 4.2 (section 4.1): Values for study-level prior s0S, s0T, rho0

# At this point we are set to fit our study-level prior. We are
# treating the observed K-M curves as the true underlying distributions
# (as opposed to noisy versions of the true underlying distributions)
# for the purposes of this simulation, so we fit the prior under this
# assumption. We will fit a prior with a zero mean. 
# mu0S=mu0T=0; optimize other 3
# Negative log-likelihood, since we will minimize
# We are optimizing in 3D: sig0S, sig0T, rho0

nLL <- function(x, type) {
  if (type == "A") {
    mu0 <- x[1:2]
  } else {
    mu0 <- c(0, 0)
  }
  sd0S <- max(1e-4, x[3])
  sd0T <- max(1e-4, x[4])
  rho0 <- min(max(x[5], -1), 1)
  Sig0 <- rbind(c(sd0S^2, rho0*sd0S*sd0T), c(rho0*sd0S*sd0T, sd0T^2))
  -sum(apply(all.sim[,c("mu.PFS", "mu.OS")], 1, function(r) {
    mui <- c(r["mu.PFS"], r["mu.OS"])
    -0.5*(as.vector((mui-mu0) %*% ginv(Sig0) %*% (mui-mu0)) + log(det(Sig0))) - log(2*pi)
  }))
}
init <- c(mean(all.sim$mu.PFS), mean(all.sim$mu.OS), sd(all.sim$mu.PFS), sd(all.sim$mu.OS), cor(all.sim$mu.PFS, all.sim$mu.OS))
optB <- optim(init, nLL, type="B")
optB$par[1:2] <- 0
print(paste("Type B params:", paste(optB$par, collapse=" ")))
print(paste("Type B AIC:", 6 + 2*optB$value, "BIC", 3*log(nrow(all.sim)) + 2*optB$value))




#######################################################################
# Simulation preliminaries: remaining study-level prior 
#   and trial design

# OK, so now we are ready to finish up with our params and actually
# run the simulation. We separately run the optimization code to fit the trial
# stopping parameters n, v1, v2, z1, z2 (see Parameterizing_prior_Simulation.R). 
# The results are hard-coded here.
params$s0S <- 0.326453619374035
params$s0T <- 0.254602972671395
params$rho0 <- 0.737315875040422
params$n <- c("A"=336, "B"=NA, "C"=451)
params$v1 = c("A"=0.02541487, "B"=NA, "C"=0.02361996)
params$v2 = c("A"=0.02135448, "B"=NA, "C"=0.02135162)
params$z1 = c("A"=3.038074, "B"=NA, "C"=2.690110)
params$z2 = c("A"=1.960735, "B"=NA, "C"=1.961836)

# Return the variance of a type C trial, as a function of its number of observations
#   nS: number of surrogate outcomes
#   nT: number of true outcomes
#   s0S: study level surrogate effect size variance
#   s0T: study level true effect size variance
#   rho0: study level correlation
#   rhoI: individual level correlation
vC <- function(nS, nT, s0S, s0T, rho0, rhoI){ 
  D0 <- s0S^2*s0T^2*(1-rho0^2)
  (s0T^2*(1-rhoI^2)+nT*D0/4+(nS-nT)*D0*(1-rhoI^2)/4) / (1-rhoI^2+nT^2*D0/16+nT*(s0T^2+s0S^2-2*rho0*rhoI*s0S*s0T)/4+(nS-nT)*(nT*D0/4+s0S^2*(1-rhoI^2))/4)
}



#########################################################################
# Functions for simulating trials using our new designs

# Let's simulate some trials. We will need to know the
# timing of our intermediate and final analyses. This differs
# by trial type (A/B/C) and target variance v.

# Return the time of analysis given:
#   time: matrix of patient event times
#   type: trial type we're analyzing
#   v: target variance
#   params: study-level parameters
analysis.time <- function(times, type, v, params) {
  if (type == "A") {
    sort(times$enroll + times$OS)[ceiling(4/v-4/params$s0T^2)]
  } else if (type == "B") {
    sort(times$enroll + times$PFS)[ceiling(4*(params$s0T^2-v)/(params$s0S^2*(v-params$s0T^2*(1-params$rho0^2))))]
  } else {
    timings <- cbind(time=c(times$enroll + times$PFS,
                             times$enroll + times$OS),
                     isOS=rep(0:1, each=nrow(times)))
    timings <- timings[order(timings[,"time"]),]
    nS <- cumsum(1-timings[,"isOS"])
    nT <- cumsum(timings[,"isOS"])
    vars <- vC(nS, nT, params$s0S, params$s0T, params$rho0, params$rhoI)
    min(timings[vars <= v,"time"])
  }
}

# Simulate an actual trial of a given type
# Returns a vector reporting different aspects of the simulated trial.
#   exp.OS.info: interpolation function for empirical distribution for OS in experiment arm
#   ctl.OS.info: interpolation function for empirical distribution for OS in control arm
#   exp.PFS.info: interpolation function for empirical distribution for PFS (or TTP) in experiment arm
#   ctl.PFS.info: interpolation function for empirical distribution for PFS (or TTP) in control arm
#   copula: copula parameters for the control and treatment arms of this study.
#   params: underlying study-level parameters
#   type: trial type we're simulating
sim.trial <- function(exp.OS.info, ctl.OS.info, exp.PFS.info,
                      ctl.PFS.info, copula, type, params) {
  # Simulate patients and compute analysis timings and hazard ratios
  for (iter in 1:100) {
    times <- sim.patients(exp.OS.info, ctl.OS.info, exp.PFS.info,
                          ctl.PFS.info, copula, params$n[type], params)
    int.time <- analysis.time(times, type, params$v1[type], params)
    if (sum(times$enroll + times$OS < int.time & times$isexp == 0) == 0 ||
        sum(times$enroll + times$OS < int.time & times$isexp == 1) == 0) {
      next
    }
    int.enroll <- sum(times$enroll <= int.time)
    int.enroll.cost <- sum(sapply(times$enroll[which(times$enroll <= int.time)], function(x){
      params$cp*exp(-params$gamma*x)
    }))
    int.ES <- effect.size(times, int.time)
    final.time <- analysis.time(times, type, params$v2[type], params)
    final.enroll <- sum(times$enroll <= final.time)
    final.enroll.cost <- sum(sapply(times$enroll[which(times$enroll <= final.time)], function(x){
      params$cp*exp(-params$gamma*x)
    }))
    final.ES <- effect.size(times, final.time)
    break
  }
  if (iter == 100) browser()  # Can't get intermed with OS events...

  # Do the Bayesian updates. In these updates, use the actual
  # variances from the Cox PH test instead of the asymptotic
  # variances.
  Sig0 <- with(params, rbind(c(s0S^2, rho0*s0S*s0T), c(rho0*s0S*s0T, s0T^2)))
  if (type == "A") {
    int.T <- rbind(c(0, 0), c(0, 1/int.ES["V.OS"]))
    int.Sig <- ginv(ginv(Sig0) + int.T)
    int.mu <- as.vector(int.Sig %*% (int.T %*% int.ES[c("mu.PFS", "mu.OS")]))
    final.T <- rbind(c(0, 0), c(0, 1/final.ES["V.OS"]))
    final.Sig <- ginv(ginv(Sig0) + final.T)
    final.mu <- as.vector(final.Sig %*% (final.T %*% final.ES[c("mu.PFS", "mu.OS")]))
  } else if (type == "B") {
    int.S <- rbind(c(1/int.ES["V.PFS"], 0), c(0, 0))
    int.Sig <- ginv(ginv(Sig0) + int.S)
    int.mu <- as.vector(int.Sig %*% (int.S %*% int.ES[c("mu.PFS", "mu.OS")]))
    final.S <- rbind(c(1/final.ES["V.PFS"], 0), c(0, 0))
    final.Sig <- ginv(ginv(Sig0) + final.S)
    final.mu <- as.vector(final.Sig %*% (final.S %*% final.ES[c("mu.PFS", "mu.OS")]))
  } else {
    int.SigI <- rbind(c(int.ES["C2.V.PFSj"], params$rhoI*sqrt(int.ES["C2.V.PFSj"]*int.ES["V.OS"])), c(params$rhoI*sqrt(int.ES["C2.V.PFSj"]*int.ES["V.OS"]), int.ES["V.OS"]))
    int.S <- rbind(c(1/int.ES["C2.V.PFSs"], 0), c(0, 0))
    int.Sig <- ginv(ginv(Sig0) + ginv(int.SigI) + int.S)
    int.mu <- as.vector(int.Sig %*% (ginv(int.SigI) %*% int.ES[c("C2.mu.PFSj", "mu.OS")] + int.S %*% c(int.ES["C2.mu.PFSs"], 0)))
    final.SigI <- rbind(c(final.ES["C2.V.PFSj"], params$rhoI*sqrt(final.ES["C2.V.PFSj"]*final.ES["V.OS"])), c(params$rhoI*sqrt(final.ES["C2.V.PFSj"]*final.ES["V.OS"]), final.ES["V.OS"]))
    final.S <- rbind(c(1/final.ES["C2.V.PFSs"], 0), c(0, 0))
    final.Sig <- ginv(ginv(Sig0) + ginv(final.SigI) + final.S)
    final.mu <- as.vector(final.Sig %*% (ginv(final.SigI) %*% final.ES[c("C2.mu.PFSj", "mu.OS")] + final.S %*% c(final.ES["C2.mu.PFSs"], 0)))
  }
  
  # Log the final decision, costs, and early stopping.
  int.v <- int.Sig[2,2]
  int.z <- int.mu[2] / sqrt(int.v-(int.v/params$s0T)^2)
  final.v <- final.Sig[2,2]
  final.z <- final.mu[2] / sqrt(final.v-(final.v/params$s0T)^2)
  if (!is.finite(int.z)) browser()
  if (abs(int.z) > params$z1[type]) {
    reject <- 1
    analysis <- 1
    enroll <- int.enroll
    time <- int.time
    enroll.cost <- int.enroll.cost
  } else {
    reject <- as.numeric(abs(final.z) > params$z2[type])
    analysis <- 2
    enroll <- final.enroll
    time <- final.time
    enroll.cost <- final.enroll.cost
  }
  
  # Return
  c(c("int.mu"=int.mu[2], "int.v"=int.v, "int.z"=int.z, "int.time"=int.time,
      "int.enroll.cost"= int.enroll.cost, "int.enroll"=int.enroll, "final.mu"=final.mu[2],
      "final.v"=final.v, "final.z"=final.z, "final.time"=final.time,
      "final.enroll.cost"= final.enroll.cost, "final.enroll"=final.enroll, "reject"=reject,
      "analysis"=analysis, "enroll"=enroll, "time"=time, "enroll.cost" = enroll.cost, "iter"=iter),
    setNames(int.ES, paste0("iES.", names(int.ES))),
    setNames(final.ES, paste0("fES.", names(final.ES))))
}




##################################################
# Execute simulations for each trial & trial type

set.seed(144)
for (rcode in names(big.res)) {
  print(rcode)
  x <- big.res[[rcode]]
  resA <- as.data.frame(t(replicate(1000, sim.trial(x$exp.OS.KM, x$ctl.OS.KM, x$exp.surr.KM, x$ctl.surr.KM, x$copula, "A", params))))
  # resB <- as.data.frame(t(replicate(1000, sim.trial(x$exp.OS.KM, x$ctl.OS.KM, x$exp.surr.KM, x$ctl.surr.KM, x$copula, "B", params))))
  resC <- as.data.frame(t(replicate(1000, sim.trial(x$exp.OS.KM, x$ctl.OS.KM, x$exp.surr.KM, x$ctl.surr.KM, x$copula, "C", params))))
  big.res[[rcode]]$resA <- resA
  # big.res[[rcode]]$resB <- resB
  big.res[[rcode]]$resC <- resC
}




########################################################################
# Analyze simulation results

# How do our simulated results compare to the large-scale version?
all.sim$zscore <- (all.sim$mu.PFS - params$s0S/params$s0T*params$rho0*all.sim$mu.OS) / (sqrt(1-params$rho0^2)*params$s0S)
# all.sim$zscore2 <- (all.sim$mu.PFS - (-0.117330337546939+0.251065684985415/0.11978555471523*0.646205714361927*all.sim$mu.OS)) / (sqrt(1-0.646205714361927^2)*0.251065684985415)
all.sim$rcode <- rownames(all.sim)
all.sim$mu.OS40 = all.sim40$mu.OS
all.sim$mu.PFS40 = all.sim40$mu.PFS
full.sim <- rbind(all.sim, all.sim)
full.sim$type.paper <- rep(c(" Type A", "Type C"), each=nrow(all.sim))
full.sim$type.pres <- rep(c(" T-Only", "Combined"), each=nrow(all.sim))
full.sim$propReject <- c(sapply(big.res, function(x) mean(x$resA$reject)),
                         sapply(big.res, function(x) mean(x$resC$reject)))



#Figure 6 in paper
pdf("../../Plots/rejectRate40.pdf", width=6, height=3)
print(ggplot(full.sim, aes(x=mu.OS40, y=propReject)) +
        facet_grid(cols=vars(type.paper)) +
        geom_point() +
        stat_smooth(method="loess", span=0.5, alpha=0.3) +
        geom_vline(xintercept=0, lty=2, col="red") +
        geom_vline(xintercept=0.5, lty=2, col="red") +
        geom_vline(xintercept=-0.5, lty=2, col="red") +
        geom_hline(yintercept=0.05, lty=2, col="red") +
        geom_hline(yintercept=0.8, lty=2, col="red") +
        xlab("True OS Hazard Ratio (20 Month Average Follow-Up)") +
        ylab("Null Hypothesis Reject Rate") +
        theme_bw())
dev.off()

  

save(big.res, file="../../Simulations/bigres.Robj")

k <- 20
all.sim$dist0 <- abs(all.sim$mu.OS40)
distLim <- sort(all.sim$dist0)[k]
nullStudies <- all.sim$rcode[all.sim$dist0 <= distLim]
tapply(full.sim$propReject[full.sim$rcode %in% nullStudies],
       full.sim$type.paper[full.sim$rcode %in% nullStudies], t.test)
all.sim$distA <- abs(abs(all.sim$mu.OS)-0.5)
distLim <- sort(all.sim$distA)[k]
altStudies <- all.sim$rcode[all.sim$distA <= distLim]
tapply(full.sim$propReject[full.sim$rcode %in% altStudies],
       full.sim$type.paper[full.sim$rcode %in% altStudies], t.test)


# Two options for prior estimation
library(mixtools)
mu <- c(mean(all.sim$mu.PFS), mean(all.sim$mu.OS))
sds <- c(sd(all.sim$mu.PFS), sd(all.sim$mu.OS))
rho <- cor(all.sim$mu.PFS, all.sim$mu.OS)
ellipse2 <- mixtools:::ellipse(mu, rbind(c(params$s0S^2, params$rho0*params$s0S*params$s0T), c(params$rho0*params$s0S*params$s0T, params$s0T^2)), draw=FALSE)
ellipse2 <- data.frame(mu.PFS=ellipse2[,1], mu.OS=ellipse2[,2])



#Figure 5 in paper
pdf("../../Plots/studyLevelwPrior.pdf", width=4, height=4)
print(ggplot(all.sim, aes(x=mu.PFS, y=mu.OS)) +
        geom_point() +
        geom_path(data=ellipse2, lty=2) +
        xlab(expression("Surrogate Effect Size" ~ (mu[S]))) +
        ylab(expression("True Outcome Effect Size" ~ (mu[T]))))
dev.off()
#check bivariate normal
HZ.test(cbind(all.sim$mu.PFS40,all.sim$mu.OS40))


#### What is the PFS / OS correlation at final analysis?
c.PFS <- sapply(big.res, function(x) cor(x$resC$fES.mu.PFS, x$resC$fES.mu.OS))
hist(c.PFS)
mean(c.PFS)
c.PFSoldj <- sapply(big.res, function(x) cor(x$resC$fES.mu.PFSj, x$resC$fES.mu.OS))
hist(c.PFSoldj)
mean(c.PFSoldj)
c.PFSj <- sapply(big.res, function(x) cor(x$resC$fES.C2.mu.PFSj, x$resC$fES.mu.OS))
hist(c.PFSj)
mean(c.PFSj)
c.bothPFS <- sapply(big.res, function(x) cor(x$resC$fES.C2.mu.PFSj, x$resC$fES.C2.mu.PFSs))
hist(c.bothPFS)
summary(c.bothPFS)
muJ <- sapply(big.res, function(x) mean(x$resC$fES.C2.mu.PFSj))
muS <- sapply(big.res, function(x) mean(x$resC$fES.C2.mu.PFSs))
cor(muJ, muS)
plot(muJ, muS)
abline(0, 1)

# Costs
full.sim$enroll <- c(sapply(big.res, function(x) mean(x$resA$enroll)),
                     sapply(big.res, function(x) mean(x$resC$enroll)))
full.sim$time <- c(sapply(big.res, function(x) mean(x$resA$time)),
                   sapply(big.res, function(x) mean(x$resC$time)))

full.sim$disc_enrollment_accurate <- c(sapply(big.res, function(x) mean(x$resA$enroll.cost)),
                                       sapply(big.res, function(x) mean(x$resC$enroll.cost)))


full.sim$disc_wait <- c(sapply(big.res, function(x) {
    time_discount <- (1-exp(-params$gamma*x$resA$time))/params$gamma
    mean(params$cw *time_discount)
  }), sapply(big.res, function(x) {
    time_discount <- (1-exp(-params$gamma*x$resC$time))/params$gamma
    mean(params$cw *time_discount)
  }))



# Information for Table 4.3

# Output: average enrollment by trial type
tapply(full.sim$enroll, full.sim$type.paper, mean)

# Output: average trial length by trial type
tapply(full.sim$time, full.sim$type.paper, mean)

# Output: undiscounted avg cost
params$cp*tapply(full.sim$enroll, full.sim$type.paper, mean) + params$cw*tapply(full.sim$time, full.sim$type.paper, mean)

# Output: discounted enrollment cost
tapply(full.sim$disc_enrollment_accurate, full.sim$type.paper, mean)

#Output: discounted waiting cost
tapply(full.sim$disc_wait, full.sim$type.paper, mean)

#Output: total discounted cost
tapply(full.sim$disc_wait, full.sim$type.paper, mean)+tapply(full.sim$disc_enrollment_accurate, full.sim$type.paper, mean)





# Explore ranges of reject rates closest to Type II trials
# (Claims in 2nd to last paragraph of section 4.3)

# minimum value of proportion rejected by trial type
full.sim %>% group_by(type.paper) %>% 
  filter(mu.OS40 >= -0.6 & mu.OS40 <= -0.4) %>% summarise(min(propReject))

# maximum value of proportion rejected by trial type
full.sim %>% group_by(type.paper) %>% 
  filter(mu.OS40 >= -0.6 & mu.OS40 <= -0.4) %>% summarise(max(propReject))





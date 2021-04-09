# Section 4.1, helper function for simulation (main.Sim.R)
# Trial design parameters mentioned in section 4.1 and used in MBC simulation
# Table 4.2: last 10 entries

library(dplyr)



#############################################################################
# Fundamental functions for Type 1 / Type 2 Analysis for Survival


# Compute alpha spend at the final analysis, given v1, z1, v2, and z2
a2 <- function(v1, z1, v2, z2, p) {
  # We pulled out the "1+" term from the integral -- this is the pnorm(z1) - pnorm(-z1)
  pnorm(z1) - pnorm(-z1) + integrate(function(x) {
    dnorm(x)*(pnorm((-z2*sqrt(1/v2-1/p$s0T^2)-x*sqrt(1/v1-1/p$s0T^2))/sqrt(1/v2-1/v1))-pnorm((z2*sqrt(1/v2-1/p$s0T^2)-x*sqrt(1/v1-1/p$s0T^2))/sqrt(1/v2-1/v1)))
  }, -z1, z1)$value
}

# Compute the z2 given v1, z1, and alpha1, along with the selected v2 and target alpha
Balpha <- function(v1, z1, alpha1, v2, p) {
  uniroot(function(z2) p$alpha-alpha1-a2(v1, z1, v2, z2, p), c(-qnorm(p$alpha/2)-1e6, 10), tol=1e-6)$root
}

# We reject 1-b2 at the final analysis, given v1, z1, v2, and z2
b2 <- function(v1, z1, v2, z2, p) {
  # We pulled out the "1+" term from the integral
  min(max(1-(pnorm(z1-p$delta*sqrt(1/v1-1/p$s0T^2)) - pnorm(-z1-p$delta*sqrt(1/v1-1/p$s0T^2)) +
               integrate(function(x) {
                 dnorm(x-p$delta*sqrt(1/v1-1/p$s0T^2))*(pnorm((-z2*sqrt(1/v2-1/p$s0T^2)-x*sqrt(1/v1-1/p$s0T^2))/sqrt(1/v2-1/v1)-p$delta*sqrt(1/v2-1/v1))-pnorm((z2*sqrt(1/v2-1/p$s0T^2)-x*sqrt(1/v1-1/p$s0T^2))/sqrt(1/v2-1/v1)-p$delta*sqrt(1/v2-1/v1)))
               }, -z1, z1)$value), 0), 1)
}

# Compute the z2 needed to reject 1-beta overall, given v1, z1, beta1, and target v2 and beta
# 1-beta = (1-b1) + (1-b2)
# b2 = 1 - b1 + beta
Bbeta <- function(v1, z1, beta1, v2, p) {
  uniroot(function(z2) 1+p$beta-beta1-b2(v1, z1, v2, z2, p), c(0, p$delta*(1/v2-1/v1)/sqrt(1/v2-1/p$s0T^2)+10), tol=1e-6)$root
}

###########
# Survival cost functions
E <- function(t, l, lE, n) {
  ret <- lE/(l*n) * ifelse(t <= n/lE, l*t + exp(-l*t) - 1, l*n/lE + exp(-l*t)-exp(-l*(t-n/lE))) 
  ret[is.nan(ret)] <- 0
  ret
}
A <- function(t, muo, lo, lE, n) (E(t, lo, lE, n) + E(t, exp(muo)*lo, lE, n)) / 2
vA <- function(nT, s0T) s0T^2/(1+nT*s0T^2/4)
vB <- function(nS, s0S, s0T, rho0){
  D0 <- s0S^2*s0T^2*(1-rho0^2)
  (s0T^2+nS*D0/4)/(1+nS*s0S^2/4)
} 
vC <- function(nS, nT, s0S, s0T, rho0, rhoI){ 
  D0 <- s0S^2*s0T^2*(1-rho0^2)
  (s0T^2*(1-rhoI^2)+nT*D0/4+(nS-nT)*D0*(1-rhoI^2)/4) / (1-rhoI^2+nT^2*D0/16+nT*(s0T^2+s0S^2-2*rho0*rhoI*s0S*s0T)/4+(nS-nT)*(nT*D0/4+s0S^2*(1-rhoI^2))/4)
}
tA <- function(muT, v, n, p) {
  if (p$s0T^2 < v) {
    0  # We immediately hit the number of people needed
  } else {
    tMaxes <- exp(seq(0, log(.Machine$double.xmax), length.out=100))
    tVars <- vA(n*A(tMaxes, muT, p$lCT, p$lE, n), p$s0T)
    if (all(tVars >= v)) {
      browser()  # Should not happen -- none of our designs hit variance
    }
    tMax <- min(tMaxes[tVars < v])
    tryCatch({
      uniroot(function(tt) {
        vA(n*A(tt, muT, p$lCT, p$lE, n), p$s0T)-v
      }, c(0, tMax), extendInt = "down", maxiter=10000)$root
    }, error = function(e) { print("A") ; browser() })
  }
}
tB <- function(muS, v, n, p) {
  if (p$s0T^2 < v) {
    0  # We immediately hit the number of people needed
  } else {
    tMaxes <- exp(seq(0, log(.Machine$double.xmax), length.out=100))
    tVars <- vB(n*A(tMaxes, muS, p$lCS, p$lE, n), p$s0S, p$s0T, p$rho0)
    if (all(tVars >= v)) {
      browser()  # Should not happen -- none of our designs hit variance
    }
    tMax <- min(tMaxes[tVars < v])
    tryCatch({
      uniroot(function(tt) {
        vB(n*A(tt, muS, p$lCS, p$lE, n), p$s0S, p$s0T, p$rho0)-v
      }, c(0, tMax), extendInt = "down", maxiter=10000)$root
    }, error = function(e) { print("B") ; browser() })
  }
}
tC <- function(muS, muT, v, n, p) {
  if (p$s0T^2 < v) {
    0  # We immediately hit the number of people needed
  } else {
    tMaxes <- exp(seq(0, log(.Machine$double.xmax), length.out=100))
    tVars <- vC(n*A(tMaxes, muS, p$lCS, p$lE, n), n*A(tMaxes, muT, p$lCT, p$lE, n), p$s0S, p$s0T, p$rho0, p$rhoI)
    if (all(tVars >= v)) {
      if (abs(muT) < 100) {
        print(paste("Warning: tC hit max t with muT=", muT))
      }
      return(10^300)
    }
    tMax <- min(tMaxes[tVars < v])
    tryCatch({
      uniroot(function(tt) {
        vC(n*A(tt, muS, p$lCS, p$lE, n), n*A(tt, muT, p$lCT, p$lE, n), p$s0S, p$s0T, p$rho0, p$rhoI)-v
      }, c(0, tMax), extendInt = "down", maxiter=10000)$root
    }, error = function(e) { print("C") ; browser() })
  }
}
pA <- function(muT, v, n, p) min(p$lE*tA(muT, v, n, p), n)
pB <- function(muS, v, n, p) min(p$lE*tB(muS, v, n, p), n)
pC <- function(muS, muT, v, n, p) min(p$lE*tC(muS, muT, v, n, p), n)

# A helper method that computes an integral of our cost functions
# over a specified range of the intermediate analysis posterior mean z-score,
# based on the assumption that we will continue to a specified target variance v.
ObjHelper <- function(n, vintermed, vfinal, zmin, zmax, tt, p) {
  if (tt == "A") {
    c((p$cw/p$gamma) * integrate(function(xx) {
        w <- sapply(xx, function(x) tA(muT=x*sqrt(vintermed-vintermed^2/p$s0T^2), v=vfinal, n=n, p))
        dnorm(xx) * (1-exp(-p$gamma*w)) #exp(-p$gamma*w) is the discount factor for value xx
      }, zmin, zmax)$value,
      (p$cp*p$lE/p$gamma) * integrate(function(xx) {
        w <- sapply(xx, function(x) tA(muT=x*sqrt(vintermed-vintermed^2/p$s0T^2), v=vfinal, n=n, p))
        dnorm(xx) * (1-exp(-p$gamma*min(w, n/p$lE)))
      }, zmin, zmax)$value)
  } else if (tt == "B") {
    #print("cost wait")
    cost_wait <- (p$cw/p$gamma) * integrate(function(xx) {
      w <- sapply(xx, function(x) {
        sX <- p$s0S/(p$rho0*p$s0T) * x*sqrt(vintermed-vintermed^2/p$s0T^2)
        tB(muS=sX, v=vfinal, n=n, p)
      })
      dnorm(xx) * (1-exp(-p$gamma*w)) 
    }, zmin, zmax)$value
    #print(cost_wait)
    
    #print("cost_people")
    cost_people <- (p$cp*p$lE/p$gamma) * integrate(function(xx) {
      w <- sapply(xx, function(x) {
        sX <- p$s0S/(p$rho0*p$s0T) * x*sqrt(vintermed-vintermed^2/p$s0T^2)
        tB(muS=sX, v=vfinal, n=n, p)
      })
      dnorm(xx) * (1-exp(-p$gamma*min(w,n/p$lE)))
    }, zmin, zmax)$value
    #print(cost_people)
    
    c(cost_wait, cost_people)
  } else if (tt == "C") {
    c((p$cw/p$gamma) * integrate(function(xx) {
      w <- sapply(xx, function(x) {
        sX <- p$rho0*p$s0S/p$s0T * x*sqrt(vintermed-vintermed^2/p$s0T^2)
        tC(muS=sX, muT=x*sqrt(vintermed-vintermed^2/p$s0T^2), v=vfinal, n=n, p)
      })
      dnorm(xx) * (1-exp(-p$gamma*w))
    }, zmin, zmax)$value,
    (p$cp*p$lE/p$gamma) * integrate(function(xx) {
      w <- sapply(xx, function(x) {
        sX <- p$rho0*p$s0S/p$s0T * x*sqrt(vintermed-vintermed^2/p$s0T^2)
        tC(muS=sX, muT=x*sqrt(vintermed-vintermed^2/p$s0T^2), v=vfinal, n=n, p)
      })
      dnorm(xx) * (1-exp(-p$gamma*min(w, n/p$lE)))
    }, zmin, zmax)$value)
  }
}

# The objective value with no intermediate analysis
ObjSimple <- function(n, v, tt, p) {
  ObjHelper(n=n, vintermed=v, vfinal=v, zmin=-Inf, zmax=Inf, tt, p)
}

# The objective value with an intermediate analysis
ObjFull <- function(n, v1, z1, v2, tt, p) {
  ObjHelper(n=n, vintermed=v1, vfinal=v1, zmin=-Inf, zmax=-z1, tt, p) +
    ObjHelper(n=n, vintermed=v1, vfinal=v2, zmin=-z1, zmax=z1, tt, p) +
    ObjHelper(n=n, vintermed=v1, vfinal=v1, zmin=z1, zmax=Inf, tt, p)
}

# What is the maximum n we need to care about? This is the minimum n such that we
# are still enrolling when we hit our target v2, even if there are no events at all
# in the treatment group. Any n larger than this will have exactly the same cost,
# since those extra patients will never be enrolled
# So we need (n/2) * Pr(observe in control arm at time T=n/lE) = "o" for target o.
# For type A trials (true outcomes only), the probability is
# 1-(1-exp(-lT*n/lE))/(lT*n/lE). So we can express this as:
# an+b = exp(-cn) for:
# a = -lT/lE
# b = 2*lT/lE*o + 1
# c = lT/lE
# For type A, we have o = 4/v1 - 4/s0T^2
# From https://www.wolframalpha.com/input/?i=solve+for+x%3A+ax%2Bb+%3D+exp%28-cx%29
# this can be solved as n = r/c-b/a, where r solves the following equation:
# r*exp(r) = c*exp(b*c/a)/a
# So we need the inverse of the function r*exp(r). This is called the Lambert W
# function (https://en.wikipedia.org/wiki/Lambert_W_function) and is efficiently
# implemented in the lamW package in R.
# For Type C, this will also be an acceptable upper bound (we essentially ignore the
# existence of surrogates). For Type B, we'll have target surrogate count
# n = 4*(s0T^2-v)/(s0S^2*(v-s0T^2*(1-rho0^2)))
# Further, we'll need to switch to rates lS instead of lT. This yields:
# a = -lS/lE
# b = 2*lS/lE*n + 1
# c = lS/lE
# install.packages("lamW")
library(lamW)
nUB <- function(v, tt, p) {
  if (tt == "B") {
    aa <- -p$lCS/p$lE
    bb <- 2*p$lCS/p$lE*4*(p$s0T^2-v)/(p$s0S^2*(v-p$s0T^2*(1-p$rho0^2))) + 1
    cc <- p$lCS/p$lE
  } else {
    aa <- -p$lCT/p$lE
    bb <- 2*p$lCT/p$lE*(4/v-4/p$s0T^2) + 1
    cc <- p$lCT/p$lE
  }
  lambertW0(cc*exp(bb*cc/aa)/aa)/cc - bb/aa
}

# Given v1, z1, and v2, what is the optimal trial size n?
# If v1 and z1 are NA, then there's no intermediate analysis and v2 is the
# final analysis timing
optN <- function(v1, z1, v2, tt, p) {
  if (tt == "B") {
    nmin <- 4*(p$s0T^2-v2)/(p$s0S^2*(v2-p$s0T^2*(1-p$rho0^2)))+1e-6
  } else {
    nmin <- 4/v2-4/p$s0T^2+1e-9
  }
  nmax <- nUB(v2, tt, p)

  # Start with a grid to get a good starting value
  ntries <- seq(nmin, nmax, length.out=20)
  thisf <- function(n) {
    if (!is.na(v1) && !is.na(z1)) {
      sum(ObjFull(n, v1, z1, v2, tt, p))
    } else {
      sum(ObjSimple(n, v2, tt, p))
    }
  }
  tryvals <- sapply(ntries, thisf)
  
  # Optimize from our starting value
  oval <- optim(ntries[which.min(tryvals)], thisf,
                method="L-BFGS-B", lower=nmin, upper=nmax)
  c(n=oval$par, cost=oval$value)
}

# Run over a grid of v1 and z1 values (NULL defaults to picking gsize points)
# For each v1/z1 pair, optimizes v2, z2, and n and returns the resulting cost.
# Also runs the design with no intermediate analysis. Returns a list with
# $grid as the grid data frame, $noint as the no-intermediate results, and
# $wins as the number of wins for the intermediate versus no-intermediate.
gridRun <- function(v1, z1, tt, p, gsize=10) {
  # Parse parameters
  zstar <- -qnorm(p$alpha/2)
  vstar <- uniroot(function(v1) pnorm(zstar-p$delta*sqrt(1/v1-1/p$s0T^2)) - pnorm(-zstar-p$delta*sqrt(1/v1-1/p$s0T^2)) - p$beta, c(1e-9, p$s0T^2), tol=1e-6)$root
  
  vcrit <- p$s0T^2*(1-p$rho0^2)  # Best achievable variance for type B

  if (tt == "B" && vcrit > vstar) {
    return(list(grid=NULL, noint=NULL, wins=NA,
                message=paste("Best possible variance", vcrit,
                              "but need variance", vstar)))
  }
  
  if (is.null(v1)) {
    # Evenly spaced effective sample sizes at intermediate analysis
    v1 <- 4/(tail(head(seq(0, 4/vstar-4/p$s0T^2, length.out=gsize+2), -1), -1)+4/p$s0T^2)
  }
  if (is.null(z1)) {
    # evenly spaced probabilities of reject at intermediate analysis
    z1 <- -qnorm(tail(head(seq(0, pnorm(-zstar), length.out=gsize+2), -1), -1))
  }

  # Design with no intermediate analysis
  ostar <- optN(NA, NA, vstar, tt, p)
  cstar <- ObjSimple(unname(ostar["n"]), vstar, tt, p)
  noint <- data.frame(v=vstar, z=zstar, n=unname(ostar["n"]),
                      cw=cstar[1], cp=cstar[2], cost=cstar[1]+cstar[2])
  
  # Run our grid
  grid <- expand.grid(v1=v1, z1=z1)
  res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(idx) {
    v1 <- grid$v1[idx]
    z1 <- grid$z1[idx]
    alpha1 <- 2*pnorm(-z1)  # Proportion of alpha spent at intermediate analysis
    beta1 <- pnorm(z1-p$delta*sqrt(1/v1-1/p$s0T^2)) - pnorm(-z1-p$delta*sqrt(1/v1-1/p$s0T^2))
    if (alpha1 >= p$alpha || beta1 <= p$beta) {
      browser()  # Numerical issues!
    }
    v2 <- uniroot(function(v2) Balpha(v1, z1, alpha1, v2, p)-Bbeta(v1, z1, beta1, v2, p), c(1e-9, v1-1e-10), tol=1e-6)$root
    if (tt == "B" && v2 <= vcrit) {
      return(NULL)  # This v1/z1 yields a sufficiently low v2 that we could never finish this trial for type B
    }
    z2 <- (Balpha(v1, z1, alpha1, v2, p)+Bbeta(v1, z1, beta1, v2, p))/2
    nopt <- optN(v1, z1, v2, tt, p)
    fullobj <- ObjFull(nopt["n"], v1, z1, v2, tt, p)
    data.frame(v1=v1, z1=z1, alpha1=alpha1, beta1=beta1, v2=v2, z2=z2, n=nopt["n"],
               cw=fullobj[1], cp=fullobj[2], cost=nopt["cost"])
  }))
  
  # Return
  list(grid=res, noint=noint, wins=sum(res$cost < noint$cost))
}

# Run only the "no intermediate analysis" designs for our three trial types
noint <- function(p) {
  # What are the target variance and z-score?
  zstar <- -qnorm(p$alpha/2)
  vstar <- uniroot(function(v1) pnorm(zstar-p$delta*sqrt(1/v1-1/p$s0T^2)) - pnorm(-zstar-p$delta*sqrt(1/v1-1/p$s0T^2)) - p$beta, c(1e-9, p$s0T^2), tol=1e-6)$root
  
  do.call(rbind, lapply(c("A", "B", "C"), function(tt) {
    if (tt == "B" && vstar <= p$s0T^2*(1-p$rho0^2)) {
      data.frame(tt="B", v=NA, z=NA, n=NA, cw=NA, cp=NA, cost=NA)
    } else {
      ostar <- optN(NA, NA, vstar, tt, p)
      cstar <- ObjSimple(unname(ostar["n"]), vstar, tt, p)
      data.frame(tt=tt, v=vstar, z=zstar, n=unname(ostar["n"]),
                 cw=cstar[1], cp=cstar[2], cost=cstar[1]+cstar[2])
    }
  }))
}



# Returns optimal designs for each trial given input 
#   p: study-level parameters 
allDesigns <- function(p){
  system.time(SurvDesignA <- gridRun(v1 = NULL, z1=NULL, "A", p))
  system.time(SurvDesignB <- gridRun(v1 = NULL, z1=NULL, "B", p))
  system.time(SurvDesignC <- gridRun(v1 = NULL, z1=NULL, "C", p))
  
  bestIntDesignA <- SurvDesignA$grid[which.min(SurvDesignA$grid$cost),]  
  bestIntDesignA <- bestIntDesignA %>%
    rename(cpInt = cp, cwInt = cw, costInt = cost, nInt = n)
  bestNoIntDesignA <- SurvDesignA$noint
  relImproveA <- (min(SurvDesignA$grid$cost) / SurvDesignA$noint$cost) 
  DesignA <- unlist(c(bestIntDesignA, bestNoIntDesignA, relImprove = relImproveA))
  
  bestIntDesignB <- SurvDesignB$grid[which.min(SurvDesignB$grid$cost),]  
  bestIntDesignB <- bestIntDesignB %>%
    rename(cpInt = cp, cwInt = cw, costInt = cost, nInt = n)
  bestNoIntDesignB <- unlist(SurvDesignB$noint)
  relImproveB <- (min(SurvDesignB$grid$cost) / SurvDesignB$noint$cost) 
  DesignB <- unlist(c(bestIntDesignB, bestNoIntDesignB, relImprove = relImproveB))
  
  bestIntDesignC <- SurvDesignC$grid[which.min(SurvDesignC$grid$cost),]  
  bestIntDesignC <- bestIntDesignC %>%
    rename(cpInt = cp, cwInt = cw, costInt = cost, nInt = n)
  bestNoIntDesignC<- unlist(SurvDesignC$noint)
  relImproveC <- (min(SurvDesignC$grid$cost) / SurvDesignC$noint$cost) 
  DesignC <- unlist(c(bestIntDesignC, bestNoIntDesignC, relImprove = relImproveC))
  
  costs_and_params <- list(DesignA = DesignA, DesignB = DesignB, DesignC = DesignC)
  costs_and_params
}




################################################################################
# Parameters for optimization

#First, input exogeneous parameters
rate <- 8  # patients/month
lE <- 8
#Maximum likelihood estimates for prior params from Hazard Ratios (see mainSim.R)
mu0S <- 0
mu0T <- 0
s0S <- 0.3264536
s0T <- 0.2546030
rho0 <- 0.7373159
#Initialize IPD correlation based on BURZY0408, and patient std based on
rhoI <- 0.688
#Initialize costs based on cost literature 
cp <- 50000
cw <- 20340000 #4938200 #3500000
#type 1/type 2 parameters
alpha <- 0.05
beta <- 0.2
delta <- 0.5

# let's get some survival designs up in here!
surr.rate <- (log(2))/7.2 #Surrogate arrival rate (events/month)
true.rate <- (log(2))/21.3 
gamma <- 0.08/12 # discount factor


p <- list(s0S=s0S, s0T=s0T, delta=delta, alpha=alpha, beta=beta, cw=cw, cp=cp,
          lCS=surr.rate, lCT=true.rate, lE=lE, rho0=rho0, rhoI=rhoI, gamma = gamma)



################################################################################
# Optimize trial design based on above parameters

#survDesign <- allDesigns(p) #This is with no intermediate analysis
system.time(ggA <- gridRun(v1=NULL, z1=NULL, "A", p, gsize=20))
system.time(ggB <- gridRun(v1=NULL, z1=NULL, "B", p, gsize=20))
system.time(ggC <- gridRun(v1=NULL, z1=NULL, "C", p, gsize=20))


#clean up results for type A trial
bestIntDesignA <- ggA$grid[which.min(ggA$grid$cost),]  
bestIntDesignA <- bestIntDesignA %>%
  rename(cpInt = cp, cwInt = cw, costInt = cost, nInt = n)
bestNoIntDesignA <- ggA$noint
relImproveA <- (min(ggA$grid$cost) / ggA$noint$cost) 
DesignA <- unlist(c(bestIntDesignA, bestNoIntDesignA, relImprove = relImproveA))

#clean up results for type B trial (could wind up infeasible)
bestIntDesignB <- ggB$grid[which.min(ggB$grid$cost),]  
bestIntDesignB <- bestIntDesignB %>%
  rename(cpInt = cp, cwInt = cw, costInt = cost, nInt = n)
bestNoIntDesignB <- unlist(ggB$noint)
relImproveB <- (min(ggB$grid$cost) / ggB$noint$cost) 
DesignB <- unlist(c(bestIntDesignB, bestNoIntDesignB, relImprove = relImproveB))

#clean up results for type C trial
bestIntDesignC <- ggC$grid[which.min(ggC$grid$cost),]  
bestIntDesignC <- bestIntDesignC %>%
  rename(cpInt = cp, cwInt = cw, costInt = cost, nInt = n)
bestNoIntDesignC<- unlist(ggC$noint)
relImproveC <- (min(ggC$grid$cost) / ggC$noint$cost) 
DesignC <- unlist(c(bestIntDesignC, bestNoIntDesignC, relImprove = relImproveC))

costs_and_params <- list(DesignA = DesignA, DesignB = DesignB, DesignC = DesignC)


# final endogenous parameters:
n <- ceiling(c(bestIntDesignA$nInt, bestIntDesignB$nInt, bestIntDesignC$nInt))
v1 <- c(bestIntDesignA$v1, bestIntDesignB$v1, bestIntDesignC$v1)
v2 <- c(bestIntDesignA$v2, bestIntDesignB$v2, bestIntDesignC$v2)
z1 <- c(bestIntDesignA$z1, bestIntDesignB$z1, bestIntDesignC$z1)
z2 <- c(bestIntDesignA$z2, bestIntDesignB$z2, bestIntDesignC$z2)


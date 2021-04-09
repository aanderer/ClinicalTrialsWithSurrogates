# Code for Appendix E
# Bandit algorithm for expected success endpoint run with different parameters
# Results in section E.2

library(mvtnorm)
#library(ggplot2)
library(MASS)
library(matrixStats)  # changed package from Rfast to matrixStats, colprods to colProds, rowprods
#library(partitions)
library(gtools)
library(rpart)
library(rpart.plot)
library(dplyr)

# Fixed information
cc <- c(0, 0)  # cS and cT 
mu0 <- c(0, 0)
Sig0 <- cbind(c(1, 0.7), c(0.7, 1))
SigI <- cbind(c(1, 0.7), c(0.7, 1))



# Grid with p_ab values
del <- 0.05
vals <- as.matrix(expand.grid(muS=seq(-5, 5, del), muT=seq(-5, 5, del)))
valsPlusCC <- t(t(vals)+cc)
p0S <- pnorm(0, cc[1]+vals[,1], 1)
p1S <- 1-p0S
p0T <- pnorm(0, cc[2]+vals[,2], 1)
p1T <- 1-p0T
p00 <- apply(valsPlusCC, 1, function(x) pmvnorm(lower=c(-Inf, -Inf), upper=c(0, 0), mean=x, sigma=SigI))
p01 <- p0S-p00
p10 <- p0T - p00
p11 <- 1-p00-p01-p10
props <- cbind(p0S, p1S, p00, p01, p10, p11)
propsS <- cbind(p0S, p1S)
propsT <- cbind(p0T, p1T)



# Simulate all our bandit algorithms on a patient population
bandit <- function(mu, N, Delta, parC, cc, SigI, mu0, Sig0, vals, props, propsS, propsT, parts) {
  # Simulate counterfactuals
  ctl <- rmvnorm(N, cc, SigI)
  trt <- rmvnorm(N, cc+mu, SigI)
  info <- data.frame(ctlS = as.numeric(ctl[,1] > 0),
                     ctlT = as.numeric(ctl[,2] > 0),
                     trtS = as.numeric(trt[,1] > 0),
                     trtT = as.numeric(trt[,2] > 0))
  info$trt <- factor(2*info$trtS + info$trtT + 1, levels=1:4)  # 1: 00; 2: 01; 3: 10; 4: 11
  
  # Work sequentially
  info$O1Full <- FALSE
  info$O2Full <- FALSE
  info$O1T <- FALSE
  info$O2T <- FALSE
  info$O1S <- FALSE
  for (i in seq_len(N)) {
    # Option 1 (Ignore pipeline), Full (aka use S and T -- Type C)
    selected <- which(info$O1Full)
    onlyS <- selected[selected > i-Delta]
    hasBoth <- selected[selected <= i-Delta]
    x <- c(sum(info$trtS[onlyS] == 0), sum(info$trtS[onlyS] == 1),
           as.vector(table(info$trt[hasBoth])))  # n0S, n0T, n00, n01, n10, n11
    post <- dmvnorm(vals, mu0, Sig0) * colProds(t(props)^x)
    post[is.na(post)] <- 0
    norm <- post / sum(post)  # grid squared normalized to sum to 1
    postMu <- colSums(vals * norm)
    postVar <- colSums(t(t(vals)-postMu)^2 * norm)
    if (postMu[2] + parC*sqrt(postVar[2]*log(i)) >= 0) {
      info$O1Full[i] <- TRUE
    }
    
    # Option 2 (Consider pipeline), Full
    # Step 1: compute the update without any consideration of pipeline
    selected <- which(info$O2Full)
    onlyS <- selected[selected > i-Delta]
    hasBoth <- selected[selected <= i-Delta]
    x <- c(sum(info$trtS[onlyS] == 0), sum(info$trtS[onlyS] == 1),
           as.vector(table(info$trt[hasBoth])))  # n0S, n0T, n00, n01, n10, n11
    post <- dmvnorm(vals, mu0, Sig0) * colProds(t(props)^x)
    post[is.na(post)] <- 0
    norm <- post / sum(post)  # grid squared normalized to sum to 1
    postMu <- colSums(vals * norm)
    postVar <- colSums(t(t(vals)-postMu)^2 * norm)
    
    if (length(onlyS) == 0) {
      pipeVar <- postVar
    } else {
      # Under our posterior mean, what is the unconditional probability of each surrogate/true pair?
      post.p00 <- pmvnorm(lower=c(-Inf, -Inf), upper=c(0, 0), mean=cc + postMu, sigma=SigI)
      post.p01 <- max(pnorm(0, cc[1]+postMu[1], 1) - post.p00, 0)
      post.p10 <- max(pnorm(0, cc[2]+postMu[2], 1) - post.p00, 0)
      post.p11 <- 1 - post.p00 - post.p01 - post.p10
      
      # We know the surrogate values for everybody in the pipeline; let's compute the expected
      # number of each surrogate/true outcome pairs once the pipeline clears:
      x.pipe <- c(0, 0, x[1]*c(post.p00, post.p01)/(post.p00+post.p01), 
                  x[2]*c(post.p10, post.p11)/(post.p10+post.p11))
      post <- dmvnorm(vals, mu0, Sig0) * colProds(t(props)^x.pipe)
      post[is.na(post)] <- 0
      norm <- post / sum(post)  # grid squared normalized to sum to 1
      pipeMu <- colSums(vals * norm)
      pipeVar <- colSums(t(t(vals)-postMu)^2 * norm)
      
      # Note -- this is a sensible estimate, but biased -- pipeMu != postMu. Still, I am struggling
      # to get exact procedures (where pipeMu == postMu) to work. Also, those will be slow (require
      # enumerating partitions of the pipeline patients). Work in progress that we may be able to leave
      # as a permanent TODO :)
      
      # Step 2: compute probability of four possible outcomes for everything in pipeline
      # post.p00 <- pmvnorm(lower=c(-Inf, -Inf), upper=c(0, 0), mean=cc + postMu, sigma=SigI)
      # post.p01 <- pnorm(0, cc[1]+postMu[1], 1) - post.p00
      # post.p10 <- pnorm(0, cc[2]+postMu[2], 1) - post.p00
      # post.p11 <- 1 - post.p00 - post.p01 - post.p10
      
      # Step 3: for all (or a random subset) of the possible sets of pipeline outcomes, compute
      # what the posterior variance would be. Use this to estimate the posterior variance once
      # the pipeline clears.
      # pp <- parts[[length(onlyS)]]
      # ppinfo <- t(apply(pp, 2, function(pipe.res) {
      # x.pipe <- c(0, 0, tail(x, 4) + pipe.res)
      # post <- dmvnorm(vals, mu0, Sig0) * colProds(t(props)^x.pipe)
      #   norm <- post / sum(post)  # grid squared normalized to sum to 1
      #   postMu <- colSums(vals * norm)
      #   postVar <- colSums(t(t(vals)-postMu)^2 * norm)
      #   c(mu=unname(postMu[2]), var=unname(postVar[2]),
      #     prob=dmultinom(pipe.res, prob=c(post.p00, post.p01, post.p10, post.p11)))
      # }))
      
      # Second step: update with expected values of event counts
      # if (postMu[2] + parC*sqrt(postVar[2]*log(i)) >= 0) {
      #   info$O1Full[i] <- TRUE
      # }
    }
    tryCatch({
      if (postMu[2] + parC*sqrt(pipeVar[2]*log(i)) >= 0) {
        info$O2Full[i] <- TRUE
      }
    }, error=function(e) browser())
    

    # Option 1, True Only
    selected <- which(info$O1T)
    hasT <- selected[selected <= i-Delta]
    x <- c(sum(info$trtT[hasT] == 0), sum(info$trtT[hasT] == 1))
    post <- dmvnorm(vals, mu0, Sig0) * colProds(t(propsT)^x)
    post[is.na(post)] <- 0
    norm <- post / sum(post)  # grid squared normalized to sum to 1
    postMu <- colSums(vals * norm)
    postVar <- colSums(t(t(vals)-postMu)^2 * norm)
    if (postMu[2] + parC*sqrt(postVar[2]*log(i)) >= 0) {
      info$O1T[i] <- TRUE
    }
    
    # Option 2, True Only
    selected <- which(info$O2T)
    onlyS <- selected[selected > i-Delta]
    hasT <- selected[selected <= i-Delta]
    x <- c(sum(info$trtT[hasT] == 0), sum(info$trtT[hasT] == 1))
    post <- dmvnorm(vals, mu0, Sig0) * colProds(t(propsT)^x)
    post[is.na(post)] <- 0
    norm <- post / sum(post)  # grid squared normalized to sum to 1
    postMu <- colSums(vals * norm)
    postVar <- colSums(t(t(vals)-postMu)^2 * norm)
    if (length(onlyS) == 0) {
      pipeVar <- postVar
    } else {
      # Under our posterior mean, what is the unconditional probability of getting a positive true outcome?
      post.1 <- pnorm(0, cc[2]+postMu[2], 1)

      # Split the pipeline according to this probability
      x.pipe <- x + length(onlyS)*c(1-post.1, post.1)
      post <- dmvnorm(vals, mu0, Sig0) * colProds(t(propsT)^x.pipe)
      post[is.na(post)] <- 0
      norm <- post / sum(post)  # grid squared normalized to sum to 1
      pipeMu <- colSums(vals * norm)
      pipeVar <- colSums(t(t(vals)-postMu)^2 * norm)
    }
    tryCatch({
      if (postMu[2] + parC*sqrt(pipeVar[2]*log(i)) >= 0) {
        info$O2T[i] <- TRUE
      }
    }, error=function(e) browser())

    
    # Option 1, Surrogate Only
    hasS <- which(info$O1S)
    x <- c(sum(info$trtS[hasS] == 0), sum(info$trtS[hasS] == 1))
    post <- dmvnorm(vals, mu0, Sig0) * colProds(t(propsS)^x)
    post[is.na(post)] <- 0
    norm <- post / sum(post)  # grid squared normalized to sum to 1
    postMu <- colSums(vals * norm)
    postVar <- colSums(t(t(vals)-postMu)^2 * norm)
    if (postMu[2] + parC*sqrt(postVar[2]*log(i)) >= 0) {
      info$O1S[i] <- TRUE
    }
  }
  
  data.frame(oracle = sum(pmax(info$ctlT, info$trtT)),
             bestArm = max(sum(info$ctlT), sum(info$trtT)),
             bO1Full = sum(info$ctlT[!info$O1Full]) + sum(info$trtT[info$O1Full]),
             bO2Full = sum(info$ctlT[!info$O2Full]) + sum(info$trtT[info$O2Full]),
             bO1T = sum(info$ctlT[!info$O1T]) + sum(info$trtT[info$O1T]),
             bO2T = sum(info$ctlT[!info$O2T]) + sum(info$trtT[info$O2T]),
             bO1S = sum(info$ctlT[!info$O1S]) + sum(info$trtT[info$O1S]),
             muS = mu[1],
             muT = mu[2])
}


#tuning parameters 0, 0.5, 1, 1.5, 2, 2.5, 3
#prior variance 0.1, 1, 10 vary together sigmaT = sigmaS
grid <- expand.grid(rhoI = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                      rho0 = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                      Delta = c(10, 25, 50, 75, 100, 125, 150),
                      sigma = c(0.1, 1, 10), parC = c(0, 0.5, 1, 1.5, 2, 2.5, 3))


#################################################################################
# Run for each set of parameters; can parallelize this
for(i in 1:21){
  
  params <- params <- grid[((i-1)*847+1):(i*847),]
  
  res <- apply(params, 1, function(x) {
    nrep <- 100
    N <- 200
    parC <- x["parC"]  # TODO: tune this
    cc <- c(0, 0)  # cS and cT
    mu0 <- c(0, 0)
    Sig0 <- unname(cbind(c(x["sigma"]^2, x["sigma"]^2*x["rho0"]), c(x["sigma"]^2*x["rho0"], x["sigma"]^2)))
    SigI <- unname(cbind(c(1, x["rhoI"]), c(x["rhoI"], 1)))
    del <- 0.1
    vals <- as.matrix(expand.grid(muS=seq(-3, 3, del), muT=seq(-3, 3, del)))
    valsPlusCC <- t(t(vals)+cc)
    p0S <- pnorm(0, cc[1]+vals[,1], 1)
    p1S <- 1-p0S
    p0T <- pnorm(0, cc[2]+vals[,2], 1)
    p1T <- 1-p0T
    p00 <- apply(valsPlusCC, 1, function(x) pmvnorm(lower=c(-Inf, -Inf), upper=c(0, 0), mean=x, sigma=SigI))
    p01 <- pmax(p0S-p00, 0)
    p10 <- pmax(p0T - p00, 0)
    p11 <- 1-p00-p01-p10
    props <- cbind(p0S, p1S, p00, p01, p10, p11)
    propsS <- cbind(p0S, p1S)
    propsT <- cbind(p0T, p1T)
    set.seed(144)
    allMu <- rmvnorm(nrep, mu0, Sig0)
    res <- do.call(rbind, apply(allMu, 1, function(mu) {
      bandit(mu, N, x["Delta"], parC, cc, SigI, mu0, Sig0, vals, props, propsS, propsT)
    }))
    list(rhoI=x["rhoI"], rho0=x["rho0"], Delta=x["Delta"], res=res)
  })
  save(res, file=paste("../Simulations/ESRes", toString(i), ".robj"))
  res.summary <- do.call(rbind, lapply(res, function(x) {
    cbind(data.frame(rhoI=x$rhoI, rho0=x$rho0, Delta=x$Delta), as.data.frame(t(as.matrix(colMeans(x$res)))))
  }))
}


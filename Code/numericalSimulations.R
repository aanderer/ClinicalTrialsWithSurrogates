# Section 3.3
# Running numerical simulations based on underlying model
# Figure 3

library(dplyr)
library(fields)
library(akima)
library(reshape)
library(viridis)


# Load functions that will help simulate trial designs
source("Type1Type2SurvivalFunc_discounting.R")


# Figure out the benefit of running trial designs with certain parameter sets
run_i <- function(i){
  full_grid <- read.csv("../Data/param_sets_comp_stats_plot_int.csv")
  
  grid<- full_grid[((i*121)-120):(i*121),2:ncol(full_grid)]
  
  for(j in 1:nrow(grid)){
    print(j)
    rho0 <- grid[j,"rho0"]
    rhoI <- grid[j,"rhoI"]
    s0S <- grid[j, "s0S"]
    s0T <- grid[j, "s0T"]
    lCT <- grid[j, "lCT"]
    cw <- grid[j, "cw"]
    p <- list(s0S=s0S, s0T=s0T, delta=0.5, alpha=0.05, beta=0.3, cw=cw, cp=1,
              lCS=0.1, lCT=lCT, lE=20, rho0=rho0, rhoI=rhoI, gamma = 0.08/12)
    print("running A")
    system.time(ggA <- gridRun(v1=NULL, z1=NULL, "A", p, gsize=20))
    print("running B")
    system.time(ggB <- gridRun(v1=NULL, z1=NULL, "B", p, gsize=20))
    print("running C")
    system.time(ggC <- gridRun(v1=NULL, z1=NULL, "C", p, gsize=20))
    print("finished trials")
    
    bestIntDesignA <- ggA$grid[which.min(ggA$grid$cost),]
   
    bestIntDesignB <- ggB$grid[which.min(ggB$grid$cost),]
    
    bestIntDesignC <- ggC$grid[which.min(ggC$grid$cost),]
    
    if(is.null(bestIntDesignA)){bestIntDesignA <- data.frame(v1=NA, z1=NA, alpha1=NA, beta1=NA, v2=NA, z2=NA, n=NA,
                                                             cw=NA, cp=NA, cost=NA)}
    if(is.null(bestIntDesignB)){bestIntDesignB <- data.frame(v1=NA, z1=NA, alpha1=NA, beta1=NA, v2=NA, z2=NA, n=NA,
                                                             cw=NA, cp=NA, cost=NA)}
    if(is.null(bestIntDesignC)){bestIntDesignC <- data.frame(v1=NA, z1=NA, alpha1=NA, beta1=NA, v2=NA, z2=NA, n=NA,
                                                             cw=NA, cp=NA, cost=NA)}
    
    costA <- bestIntDesignA$cost
    costB <- bestIntDesignB$cost
    costC <- bestIntDesignC$cost
    
    
    
    if(is.na(costB) && is.na(costA)){
      best_competitor <- NULL
      relcost <- NULL
    }else{
      div <- min(costA, costB, na.rm=TRUE)
      best_competitor <- if(costA == div){"A"} else{"B"}
      relcost <- costC/div
    }
    
    
    results_row <- cbind(rho0, rhoI, s0S, s0T, lCT, cw, bestIntDesignA, bestIntDesignB, bestIntDesignC, best_competitor, relcost)
    grid[j,]<- results_row
  }
  
  write.csv(grid, paste("../Simulations/comp_stats_params_int_", toString(i), ".csv", sep=""))
}





#####################################################################################
# Run simulation for all different parameter sets (can be run in parallel)
for(i in 1:81){
  print(i)
  run_i(i)
}






####################################################################################3
# Combine results from simulation runs
finaldf <- read.csv("../Simulations/comp_stats_params_int_1.csv")

for(i in c(2:81)){
  temp <- read.csv(paste("../Simulations/comp_stats_params_int_", toString(i), ".csv", sep=""))
  finaldf <- rbind(finaldf, temp)
  
}

write.csv(finaldf, "../Simulations/comp_stats_params_int_full.csv")
finaldf_int <- read.csv("../Simulations/comp_stats_params_int_full.csv")






################################################################################
# Plot figure 3

finaldf_int$adj_rel_cost <- sapply(finaldf_int$rel_cost, function(x){
  if(x >= 0.79){return(x)}else{return(0.79)}
})

avg_relcost_int <- finaldf_int %>% filter(s0S >= s0T) %>% group_by(rho0, rhoI) %>% summarise(relcost = mean(adj_rel_cost))
avg_relcost_int$relcost <- (1-avg_relcost_int$relcost)

im <- with(avg_relcost_int,interp(rho0,rhoI,relcost))

tmp <- data.frame(im[[3]])
tmp$rho0 <- im[[1]]
names(tmp) <- c(im[[2]],"rho0")
tmp2 <- melt(tmp, id="rho0")
tmp2$variable <- as.numeric(as.character(tmp2$variable))

pdf("../Plots/comp_stats_params_int_full.pdf", width=4.5, height=3.5)
print(ggplot(tmp2, aes(rho0, variable, fill= value)) + geom_tile() + labs(fill='Value') +
        scale_fill_viridis(discrete=FALSE) +
        xlab(expression("Study-Level Correlation" ~ (rho[0]))) +
        ylab(expression("Individual-Level Correlation" ~ (rho[I]))) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.line = element_line(colour = "black")))
dev.off()



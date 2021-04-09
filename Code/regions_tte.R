# Section 3.2
# Analyzing the regions where comparative statics of type C trial behave differently
# Figure 2

rm(list=ls())
library(ggplot2)
library(reshape2)

# fixed params for this example
s0S = 1
s0T = 1.5
oS = 200
oT = 150

# intermediate params
c = oT/oS
eps = 1/oS

# Lagrange objective for Trial Type C
objC <- function(rho0, rhoI){
  Vc = eps*(4 + 4*(-1+c)*rhoI^2)/c +
    16*eps^2*( (s0S + (c-1)*rhoI^2*s0S)^2 -2*c*rho0*rhoI*(1 + (c-1)*rhoI^2)*s0S*s0T +c^2*rhoI^2*s0T^2) / (c^2*(rho0^2-1)*s0S^2*s0T^2)
  return(-Vc)
}

vals_rho0 = 1:99/100
vals_rhoI = 1:99/100



root_drho0_1 <- function(rhoI){
  val = (s0S + (-1 + c)*rhoI^2*s0S)/(c*rhoI*s0T)
  return(val)
}

root_drho0_2 <- function(rhoI){
  val = (c*rhoI*s0T)/(s0S + (-1 + c)*rhoI^2*s0S)
  return(val)
}

# d(-Vc)/drhoI=0, solve for rho0(rhoI), 1st & 2nd solns

root_drhoI_1 <- function(rhoI){
  val = (1/((-1 + c)*c*rhoI*s0S*s0T))*(2*c*eps - 6*c*eps*rhoI^2 + 6*c^2*eps*rhoI^2 + sqrt(c*(4*c*(eps + 3*(-1 + c)*eps*rhoI^2)^2 - (-1 + c)*rhoI^2*(8*(-1 + c)*eps*(1 + (-1 + c)*rhoI^2)*s0S^2 + c*(4*c*eps + s0S^2 - c*s0S^2)*s0T^2))))
  return(val)
}

root_drhoI_2 <- function(rhoI){
  val = (1/((-1 + c)*c*rhoI*s0S*s0T))*(2*c*eps - 6*c*eps*rhoI^2 + 6*c^2*eps*rhoI^2 - sqrt(c*(4*c*(eps + 3*(-1 + c)*eps*rhoI^2)^2 - (-1 + c)*rhoI^2*(8*(-1 + c)*eps*(1 + (-1 + c)*rhoI^2)*s0S^2 + c*(4*c*eps + s0S^2 - c*s0S^2)*s0T^2))))
  return(val)
}

sols <- data.frame(rho = vals_rho0, sol1 = NA, sol2 = NA, sol3 = NA, sol4 = NA)
for(i in 1:nrow(sols)){
  rho0 = sols[i,1]
  tmp <- try(uniroot(function(x) root_drho0_1(x) -rho0, lower = 0, upper = 1, tol = 1e-10)$root,
             silent = TRUE)
  if(class(tmp) == "try-error") tmp = NA
  sols[i,2] <- tmp
  tmp <- try(uniroot(function(x) root_drho0_2(x) -rho0, lower = 0, upper = 1, tol = 1e-10)$root,
             silent = TRUE)
  if(class(tmp) == "try-error") tmp = NA
  sols[i,3] <- tmp
  tmp <- try(uniroot(function(x) root_drhoI_1(x) -rho0, lower = 1e-5, upper = 1, tol = 1e-10)$root,
             silent = TRUE)
  if(class(tmp) == "try-error") tmp = NA
  sols[i,4] <- tmp
  tmp <- try(uniroot(function(x) root_drhoI_2(x) -rho0, lower = 1e-5, upper = 1, tol = 1e-10)$root,
             silent = TRUE)
  if(class(tmp) == "try-error") tmp = NA
  sols[i,5] <- tmp
}

names(sols) <- c("rho0", "drho0_1", "drho0_2", "drhoI_1", "drhoI_2")
df <- melt(sols, id="rho0")

pdf("../Plots/region1_tte.pdf", width=4, height=3.5)
print(ggplot(data=df, aes(x=rho0,y=value,group=variable)) + 
  geom_rect(mapping=aes(xmin=0, xmax=1, ymin=0, ymax=1), colour="black", fill = NA,alpha=0.5) +
  geom_line(size = 1.25, linetype = "dashed") +
    xlab(expression("Study-Level Correlation" ~ (rho[0]))) +
    ylab(expression("Individual-Level Correlation" ~ (rho[I]))) + 
  theme_bw() + theme(legend.position = "none") +
  geom_text(label="Region 1", x=0.4, y=0.7) +
  geom_text(label="Region 4", x=0.89, y=0.925) +
  geom_text(label="Region 2", x=0.65, y=0.3) +
  geom_text(label="Region 3", x=0.875, y=0.05) 
)
dev.off()

### now take s0S > s0T

# fixed params
s0S=1
s0T=1

sols <- data.frame(rho = vals_rho0, sol1 = NA, sol2 = NA, sol3 = NA, sol4 = NA)
for(i in 1:nrow(sols)){
  rho0 = sols[i,1]
  tmp <- try(uniroot(function(x) root_drho0_1(x) -rho0, lower = 0, upper = 1, tol = 1e-10)$root,
             silent = TRUE)
  if(class(tmp) == "try-error") tmp = NA
  sols[i,2] <- tmp
  tmp <- try(uniroot(function(x) root_drho0_2(x) -rho0, lower = 0, upper = 1, tol = 1e-10)$root,
             silent = TRUE)
  if(class(tmp) == "try-error") tmp = NA
  sols[i,3] <- tmp
  tmp <- try(uniroot(function(x) root_drhoI_1(x) -rho0, lower = 1e-5, upper = 1, tol = 1e-10)$root,
             silent = TRUE)
  if(class(tmp) == "try-error") tmp = NA
  sols[i,4] <- tmp
  tmp <- try(uniroot(function(x) root_drhoI_2(x) -rho0, lower = 1e-5, upper = 1, tol = 1e-10)$root,
             silent = TRUE)
  if(class(tmp) == "try-error") tmp = NA
  sols[i,5] <- tmp
}

##############################################
# Figure 2

names(sols) <- c("rho0", "drho0_1", "drho0_2", "drhoI_1", "drhoI_2")
df <- melt(sols, id="rho0")
pdf("../Plots/region2_tte.pdf", width=4, height=3.5)
print(ggplot(data=df, aes(x=rho0,y=value,group=variable)) + 
        geom_rect(mapping=aes(xmin=0, xmax=1, ymin=0, ymax=1), colour="black", fill = NA,alpha=0.5) +
        geom_line(size = 1.25, linetype = "dashed") +
        xlab(expression("Study-Level Correlation" ~ (rho[0]))) +
        ylab(expression("Individual-Level Correlation" ~ (rho[I]))) + 
        theme_bw() + theme(legend.position = "none") +
        geom_text(label="Region 1", x=0.25, y=0.75) +
        geom_text(label="Region 2", x=0.6, y=0.4) +
        geom_text(label="Region 3", x=0.85, y=0.075) )
dev.off()

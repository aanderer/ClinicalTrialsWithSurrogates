# Section 3.3
# Look at study-level vs. individual level correlation for different surrogate/true
#   outcome pairs
# Figure 4

library(ggplot2)
library(dplyr)


########################################3
# Load data
ipd_correlations <- read.csv("../Data/IPD_surrogate_correlations.csv")

ipd_correlations <- ipd_correlations %>% filter(!is.na(Treatment.Effect.Correlation..Rho_0.) & !is.na(Numerical.Relationship.to.OS..Rho_I.))


###################################################
# Plot figure 4
png("../Plots/actual_rhoI_vs_rho0.png")
ggplot(ipd_correlations, aes(x = abs(Treatment.Effect.Correlation..Rho_0.), y = abs(Numerical.Relationship.to.OS..Rho_I.))) + geom_point() + 
  xlab(expression("Study-Level Correlation" ~ (rho[0]))) +
  ylab(expression("Individual-Level Correlation" ~ (rho[I]))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dev.off()



####################################################
# Which diseases are in our areas of interest?

# low rho0, high rhoI
ipd_correlations %>% filter(Treatment.Effect.Correlation..Rho_0. < 0.3 & Numerical.Relationship.to.OS..Rho_I. > 0.6) %>%
  select(Paper, Disease, Surrogate, End.Outcome)

# moderate rho0, low rhoI
ipd_correlations %>% filter(Treatment.Effect.Correlation..Rho_0. < 0.7 & Treatment.Effect.Correlation..Rho_0. > 0.4 & Numerical.Relationship.to.OS..Rho_I. < 0.5) %>%
  select(Paper, Disease, Surrogate, End.Outcome)


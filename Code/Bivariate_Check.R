#Code for Appendix F.1 and corresponding Figure 7
library(stringr)
library("mvnTest")
library(ggplot2)
library(dplyr)

#read in data
dat <- read.csv("../data/Bivariate Normal Check.csv", stringsAsFactors = FALSE)
dat$Paper <- str_replace_all(dat$Paper, "[\n]", "")
dat$logHR.Surrogate <- log(dat$HR.Surrogate)
dat$logHR.True <- log(dat$HR.True)

# HZ test for each paper/surrogate pair (true outcome is OS HR in all cases)
HZ <- as.data.frame(dat %>%
  group_by(Paper, Surrogate) %>%
  summarize(True=first(TRUE.), pval=HZ.test(cbind(logHR.Surrogate, logHR.True))@p.value))
print(paste("Pairs", nrow(HZ), "studies", length(unique(HZ$Paper)),
            "reject", sum(HZ$pval < 0.05), "don't reject", sum(HZ$pval >= 0.05)))

# Store pval rank of Paper/Surrogate pair in dat
pvalrank <- setNames(rank(HZ$pval), paste(HZ$Paper, HZ$Surrogate))
dat$rank <- pvalrank[paste(dat$Paper, dat$Surrogate)]
last.accept <- sum(HZ$pval < 0.05)
dat$color <- sapply(dat$rank, function(x){if(x <= last.accept){"Reject"}else{"Don't Reject"}})

# Plot effect sizes, ordered by pval
pdf("../Plots/bivariate_check.pdf", width = 6.5, height=9)
print(ggplot(dat, aes(x=logHR.Surrogate, y=logHR.True, color = color))+
  geom_point() +
  facet_wrap(vars(rank), nrow=6, ncol=5) +
  xlab("Surrogate Effect Size (Outcome Varies by Pane)") +
  ylab("ln(Overall Survival Hazard Ratio)") +
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position="none") +
  scale_color_manual(breaks = c("Reject", "Don't Reject"),
                     values=c("black", "red")))
dev.off()

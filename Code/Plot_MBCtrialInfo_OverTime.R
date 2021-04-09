# Appendix F.2
# Checking Equipoise assumption
# Figure 8

library(stringr)
library(dplyr)

OS_HR <- read.csv('../Data/OS_HR.csv')


OS_HR$year <- sapply(OS_HR$rcode, function(x){
  x <- str_trim(toString(x))
  end <- nchar(x)
  if(substr(x,end,end) == "A" | substr(x,end,end) == "B" | substr(x,end,end) == "C"){
    x <- substr(x, 1, end-2)
    end <- nchar(x)}
  year <- as.integer(substr(x, end-1, end))
  if(is.na(year)){print(x)}
  else if(year < 50){
    year <- 2000+year
  } else{year <- 1900+year}
  year
})

OS_HR$month <- sapply(OS_HR$rcode, function(x){
  x <- str_trim(toString(x))
  end <- nchar(x)
  if(substr(x,end,end) == "A" | substr(x,end,end) == "B" | substr(x,end,end) == "C"){
    x <- substr(x, 1, end-2)
    end <- nchar(x)}
  month <- as.integer(substr(x, end-3, end-2))
})

OS_HR$year[341] <- 1984
OS_HR$month[341] <- 8
OS_HR$year[422] <- 2010
OS_HR$month[422] <- 3




plot(OS_HR$year, OS_HR$effect)


#Figure 8
png("../Plots/effectSizesOverTime.png")
ggplot(OS_HR, aes(x = year, y = effect)) + geom_point() + stat_smooth() + xlab("Year of Study") + ylab("OS Effect size (log HR)")
dev.off()


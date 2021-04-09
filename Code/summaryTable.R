# Code for Section 4.1, 
# Aggregate patient characteristics for trials included in our simulation
# Table 4.1

# hard-coded rcodes for all papers used in simulation. See mainSim.R (line 111)
keep.rcode <- c("ACKLA0201 ",  "ALBAI0808 ",  "ANDER0199 ",  "ARPIN0203 ",  "BACHE0711 ", 
                "BACHE0812 ",  "BARRI0510 ",  "BASEL0512 ",  "BATIS0301 ",  "BENNE1088 ", 
                "BERGH0312 ",  "BERGH0612 ",  "BISHO0899 ",  "BLACK0310 ",  "BLAJM0399 ", 
                "BLOHM0710 ",  "BLOMQ0393 ",  "BOCCA1285 ",  "BONNE1004 ",  "BONNE1102 ", 
                "BONTE0698 ",  "BONTE1005 ",  "BRUFS0811 ",  "BUZDA0696 A", "BUZDA0696 B",
                "BYRNE0997 ",  "CHAN0409 ",   "CONTE0387 ",  "CONTE0796 ",  "CONTE0804 ", 
                "DELEN0290 ",  "DELMA0401 ",  "DICKL0816 ",  "DIERA0495 ",  "EJLER0604 ", 
                "FRENC0488 ",  "GENNA0806 ",  "GOSS0199 ",   "GRADI0113 ",  "HAMBE0411 ", 
                "HARBE0117 ",  "HARBE0316 ",  "HARRI0102 ",  "HATSC0212 ",  "HENDE0589 ", 
                "ICLI0205 ",   "INOUE0110 ",  "JASSL0301 ",  "JOHNS1109 ",  "KAUFM0400 ", 
                "KAUFM1109 ",  "LANG0213 ",   "LANGR1105 ",  "LEO1010 ",    "MAVRO0110 ", 
                "MEHTA0812 ",  "MILES0117 ",  "MILLN1207 ",  "MILLO0205 ",  "MUSS0385 ",  
                "NAMER0601 ",  "NIELS0600 ",  "NIELS1211 ",  "OBRIE0304 ",  "OSHAU0111 ", 
                "OSHAU0602 ",  "OSHAU0901 ",  "PACIL0506 ",  "PAGAN0217 ",  "PALLI0512 ", 
                "PAPAD0309 ",  "PARIT0200 ",  "PARK0513 ",   "PARK1117 ",   "PAVES0395 ", 
                "POWLE0684 ",  "RIVER0408 ",  "ROBER0411 ",  "ROBER0606 ",  "SJOST0899 ", 
                "SLAMO0301 A", "SLAMO0301 B", "SPARA0710 ",  "SPARA0909 ",  "SPARA1204 ", 
                "STEMM0310 ",  "STEMM0310b ", "VALER0111 ",  "VERMA1112 ",  "WILDI0910 ", 
                "WOLFF0113 ",  "YAMAM0217 ",  "ZHANG0317 ")
dat <- read.csv("../Data/MBCinfo.csv")
trim.trailing <- function (x) sub("\\s+$", "", x)
dat$rcode <- paste(trim.trailing(as.character(sapply(strsplit(as.character(dat$Unique_ID), "_"), "[", 1))), dat$Rand_Group)
dat$meanECOG <- 0*dat$ECOG_0 + 1*dat$ECOG_1 + 2*dat$ECOG_2 + 3*dat$ECOG_3 + 4*dat$ECOG_4

# Summary of publication years of simulated studies (row 2 table 4.1)
spl <- split(dat[dat$rcode %in% keep.rcode,], dat$rcode[dat$rcode %in% keep.rcode])
pubYear <- sapply(spl, function(x) sum(x$N_Patient*x$Pub_Year) / sum(x$N_Patient))
summary(pubYear)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1984    2000    2006    2005    2011    2017 

# Summary of number of patients in simulated studies (row 3 table 4.1)
nPatient <- sapply(spl, function(x) sum(x$N_Patient))
summary(nPatient)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.0   163.0   236.0   319.7   405.0  1286.0 

# Summary of proportion female in simulated studies
fracMale <- sapply(spl, function(x) sum(x$N_Patient*x$FRAC_MALE) / sum(x$N_Patient))
summary(fracMale)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# 0.000000 0.000000 0.000000 0.000160 0.000000 0.008947        9 
summary(1-fracMale)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#   0.9911   1.0000   1.0000   0.9998   1.0000   1.0000        9 


# Summary of median age in simulated studies
medianAge <- sapply(spl, function(x) sum(x$N_Patient*x$AGE_MED) / sum(x$N_Patient))
summary(medianAge)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   47.50   53.50   55.49   56.46   58.02   71.68       3 

# Summary of mean ECOG performance status in simulated studies
meanECOG <- sapply(spl, function(x) sum(x$N_Patient*x$meanECOG) / sum(x$N_Patient))
summary(meanECOG)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.2088  0.4330  0.5575  0.5889  0.7124  1.1248      29 

# Summary of proportion with visceral disease in simulated studies
propVisc <- sapply(spl, function(x) sum(x$N_Patient*x$Prop_Visceral) / sum(x$N_Patient))
summary(propVisc)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.1839  0.5292  0.6602  0.6312  0.7330  0.8741      34 

# Summary of reported median Overall Survival
MOS <- sapply(spl, function(x) sum(x$N_Patient*x$OS) / sum(x$N_Patient))
summary(MOS)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   9.016  17.061  20.149  21.893  26.289  45.922       9

#  Summary of reported Overall Survival log Hazard Ratios
OS.HR <- subset(read.csv("../Data/OS_HR.csv", stringsAsFactors = FALSE), rcode %in% keep.rcode)
OS.HR$exptype <- dat$Arm_Type[match(OS.HR$IDexp, dat$Unique_ID)]
OS.HR$ctltype <- dat$Arm_Type[match(OS.HR$IDctl, dat$Unique_ID)]
OS.HR$true <- ifelse(OS.HR$exptype == "experiment", OS.HR$effect, -OS.HR$effect)
summary(OS.HR$true)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.79851 -0.20333 -0.10536 -0.05225  0.07519  1.00063 

# Summary of reported median Progression Free Survival
MPFS <- sapply(spl, function(x) sum(x$N_Patient*x$PFS) / sum(x$N_Patient))
summary(MPFS)
#    Min.   1st Qu.   Median    Mean  3rd Qu.    Max.    NA's 
#   2.240     5.550    7.677   7.694    9.423  17.920      47 

# Summary of reported Progression free Survival Hazard Ratios
PFS.HR <- subset(read.csv("../Data/PFS_HR.csv", stringsAsFactors = FALSE), rcode %in% keep.rcode)
TTP.HR <- subset(read.csv("../Data/TTP_HR.csv", stringsAsFactors = FALSE), rcode %in% keep.rcode & !(rcode %in% PFS.HR$rcode))
PFS.HR <- rbind(PFS.HR, TTP.HR)
PFS.HR$exptype <- dat$Arm_Type[match(PFS.HR$IDexp, dat$Unique_ID)]
PFS.HR$ctltype <- dat$Arm_Type[match(PFS.HR$IDctl, dat$Unique_ID)]
PFS.HR$true <- ifelse(PFS.HR$exptype == "experiment", PFS.HR$effect, -PFS.HR$effect)
summary(PFS.HR$true)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.85373 -0.39244 -0.16142 -0.14234  0.04583  1.23906 

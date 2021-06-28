# ClinicalTrialsWithSurrogates

This repository contains data on Surrogate outcomes in clinical trials. The data is separated into two categories. 
The first category is data on different surrogate and true outcome pairs considered in oncology.
The second category is data specifically from Metastatic Breast Cancer trials, including data on both true and surrogate outcomes. 
The data is described in more detail below.


## Citing this work
The code and data included here are used for data analysis and simulation in the paper "Adaptive Clinical Trial Design with Surrogates: When Should We Bother?"

If you use the data or code included here, please use the following citation:



```
@article{anderer2021adaptive,
  title={Adaptive Clinical Trial Designs with Surrogates: When Should We Bother?},
  author={Anderer, Arielle and Bastani, Hamsa and Silberholz, John},
  journal={Management Science},
  year={2021}
}
```

## Data

### Different time-to-event surrogates
This group of data files contains data from meta-analyses on time-to-event surrogate outcomes, mostly in oncology trials. An example of how this data is useful can be found in the following plot:

<img align="center" width="400" src="Plots/actual_rhoI_vs_rho0.png">

Here, each dot represents a different surrogate/true outcome pair. We can see the individual-level correlation versus the study-level correlation for each pair. This helps us visualize what the relationship between these two values looks like across different potentially useful surrogate/true outcome pairs.


<b>TRUE_SURR_HR_Distribution.csv</b>: This file has data at the trial level on 30 different surrogate/true outcome pairs (we have information on between 7 and 36 trials for each pair) about hazard ratios for both outcomes. This information was compiled from scraping data from 21 meta-analyses for different diseases with time-to-event outcomes.

* Paper: id for meta-analysis it comes from (first 4 letters of first author’s last name + mmyy of publication)
* Disease: disease this trial is looking to treat
* TRUE: name or abbreviation of true outcome measured by the trial
* Surrogate: name or abbreviation of surrogate outcome measured by the trial
* Trial: Trial name/id
* Patients Treatment: treatment group size
* Patients Control: control group size
* HR Surrogate: the hazard ratio measured for the surrogate outcome
* CI Surrogate: confidence interval for HR Surrogate
* HR True: the hazard ratio measured for the true outcome 
* CI True: confidence interval for HR True


IPD_surrogate_correlations.csv: This file has information at the meta-analysis level on different surrogate/true outcome pairs about individual level and study level correlations. Compiled from scraping data from meta-analyses for different diseases with time-to-event outcomes. 
To collect these papers, we searched PubMed for meta-analyses of surrogate time-to-event endpoints and wound up with 80 papers. 
Some of these papers contain multiple surrogate endpoints for the same disease.

* Paper: id for meta-analysis it comes from (first 4 letters of first author’s last name + mmyy of publication)
* Disease: disease this trial is looking to treat
* Surrogate: name or abbreviation of surrogate outcome measured in the meta-analysis
* End Outcome: name or abbreviation of true outcome of interest measured in the meta-analysis
* Numerical Relationship to OS (Rho_I): reported individual level correlation (if
included)
* Rho_I 95% CI: reported confidence interval for individual level correlation (if
included)
* Treatment Effect Correlation (Rho_0): reported study level correlation (if
included)
* Rho_0 95% CI: reported confidence interval for study level correlation (if
included)
* Notes: any additional concerns



### Metastatic Breast Cancer (MBC)

This is historical trial-specific data for different MBC drug therapies. All of this data is pulled manually from repository of 1,865 studies of MBC drug therapies collected by Silberholz et al. (2019), which is publicly available at http://www.cancertrials.info


KM curves: This folder contains Kaplan-Meier curves pulled manually from the papers in the Silberholz et al. (2019) repository. It contains both images and data extracted from these images.

* .png files: one for each Kaplan Meier curve found in the papers in the
MBC repository. Labeled as “paper id”_”outcome”.png
* .xml files: one for each arm in every Kaplan Meier curve .png file.
Labeled as “paper id”_”outcome”_”arm id”.xml

MBC_data_final.csv: This file has information collected for each arm of each trial in MBC repository on the number of patients, the reported Hazard Ratios and related information for the different endpoints, and the names of the Kaplan-Meier curve data files if they exist.







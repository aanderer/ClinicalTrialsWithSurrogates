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

Here, each dot represents a different surrogate/true outcome pair. We can see the individual-level correlation.


<b>Bivariate Normal Check.csv</b>: This file has information from 30 different surrogate/true outcome pairs at the trial level (we have information on between 7 and 36 trials for each pair) about hazard ratios for both outcomes. This information was compiled from scraping data from 21 meta-analyses for different diseases with time-to-event outcomes.

* Paper: id for meta-analysis it comes from (first 4 letters of first author’s last name + mmyy of publication)
* Disease: disease this trial is looking to treat
* TRUE: which true outcome measured by the trial
* Surrogate: which surrogate outcome measured by the trial
* Trial: Trial name/id
* Patients Treatment: treatment group size
* Patients Control: control group size
* HR Surrogate: the hazard ratio measured for the surrogate outcome
* CI Surrogate: confidence interval for HR Surrogate
* HR True: the hazard ratio measured for the true outcome 
* CI True: confidence interval for HR True


IPD_surrogate_correlations.csv: This file has information on different surrogate/true outcome pairs. Compiled from scraping data from meta-analyses for different diseases with time-to-event outcomes. 
To collect these papers, we searched PubMed for meta-analyses of surrogate time-to-event endpoints and wound up with 80 papers. 
Some of these papers contain multiple surrogate endpoints for the same disease.



### Metastatic Breast Cancer








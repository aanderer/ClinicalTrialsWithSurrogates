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

### Different potential surrogate and true outcome pairs in oncology trials
This section 

<img align="center" width="600" src="Plots/actual_rhoI_vs_rho0.png">


Bivariate Normal Check.csv: Compiled from scraping data from 30 meta-analyses for different diseases with time-to-event outcomes.

* Paper: id for meta-analysis it comes from (first 4 letters of first author’s last name + mmyy of publication)
* Disease: disease this trial is looking to treat


IPD_surrogate_correlations.csv: Compiled from scraping data from meta-analyses for different diseases with time-to-event outcomes. 
To collect these papers, we searched PubMed for meta-analyses of surrogate time-to-event endpoints and wound up with 80 papers. 
Some of these papers contain multiple surrogate endpoints for the same disease.



### Metastatic Breast Cancer








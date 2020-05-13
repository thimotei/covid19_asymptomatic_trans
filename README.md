## The contribution of asymptomatic SARS-CoV-2 infections to transmission - a model-based analysis of the Diamond Princess outbreak

This repository contains all the data and files necessary to re-run the analyses in:

_Emery JC, Russel TW, Liu Y, et al. The contribution of asymptomatic SARS-CoV-2 infections to transmission - a model-based analysis of the Diamond Princess outbreak. medRxiv Published Online First: 11 May 2020. doi:10.1101/2020.05.07.20093849_


Full details of the study can be found in the main paper and supplementary materials, either on medRxiv or the CMMID repository:

medRxiv: https://www.medrxiv.org/content/10.1101/2020.05.07.20093849v1

CMMID repository: https://cmmid.github.io/topics/covid19/asymp-transmission.html

### Repository structure 

Contents of the folders in the repository:

`data`: Raw and processed data 

`bi`: LibBi model files

`R`: Files for processing data, executing MCMC runs, sampling the resultant posteriors to produce model outputs and plotting/exporting to tables

`posteriors`: posteriors produced by MCMC runs 

### Implementation 

To perform analyses as per the main paper or supplementary materials:

1. Raw data is stored in the folder `data`, which is then processed using `data_prep.R`. The processed files are saved back to the folder `data`. 

2. A MCMC run is then executed with `run_libbi.R` from the folder `R`, which uses the processed data and a LibBi model file from the folder `bi`. The LibBi model file used for the primary analysis is `primary.bi`. The resultant posterior file is then saved to the folder `posteriors`, where the posterior for the primary analysis is `primary.rds`.  

3. The resultant posterior is then sampled using `model_outputs.R` from the folder `R`, which requires the appropriate posterior file and the corresponding 'sampler' version of the LibBi model file used in the MCMC run. For the primary analysis these are `primary.rds` and `primary_sampler.bi` respectively. After sampling, model outputs are then calculated and stored locally in the R session.

4. `plots_and_tables.R` is then run in the same session, using the model outputs to produce the figures and tables found in both the main paper and supplementaty materials. 

The full primary analysis can be re-run by following steps 1-4 in the above using the files as they are currently set up in the repository.

The results from the primary analysis can be re-processed and re-plotted without having to peform an MCMC run by following steps 3-4 in the above and using the files as they are currently set up in the repository.

### Further resources 

LibBi is used to perform the MCMC runs and sample the resultant posteriors. The R packages RBi and RBi.helpers are used to interact with LibBi using R.

RBi documentation: https://cran.r-project.org/web/packages/rbi/vignettes/introduction.html

RBi.helpers documentation: http://sbfnk.github.io/rbi/rbi_helpers.html

LibBi documentation: https://libbi.org/documentation.html

 









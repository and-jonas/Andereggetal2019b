# Disease phenotyping for resistance breeding

This repository contains functions and analysis scripts required to reproduce results reported by Anderegg *et al*. (2019). 

## Author


> Jonas Anderegg  
> Crop Science Group  
> ETH Zürich  


## Dependencies
Scripts build primarily on functions from the R packages `prospectr`, `caret` and `tidyverse`.

## Content  

Folder `Utils` contains functions for the pre-processing of spectra, calculation and evaluation of spectral indices, training and evaluation of full-spectrum models, recursive feature elimination and analysis of variance. 

Folder `Analysis`contains scripts to implement the analysis and obtain results contained in the study. 

### Utils

#### 001_spectra_utils.R

1. `f_spc_smooth` applies the Savitzky-Golay smoothing filter to raw spectra (wrapper for `prospectr::SavitzkyGolay`)
2. `f_spc_avg` averages replicate measurements
3. `f_calc_si`calculates spectral indices
4. `f_scale_si` scales spectral indices and predictions obtained from full-spectrum model to range from 0 to 10
5. `f_spc_bin` computes average values of a signal in pre-determined bins (wrapper for `prospectr::binning`)
6. `f_spc_trim` removes noisy parts of the signal in pre-determined ranges
7. Several data wrangling helper functions. 

#### 002_pls_utils.R

1. `f_splsda` fits PLS-DA and S-PLS-DA models for classification into healthy and disease canopies. This is essentially a wrapper for functions of the `mixOmics` package. 
2. `assess_mod_perf` evaluates model performance on the main experiment.
3. `assess_mod_perf2` evaluates model performance on FPWW023 (used to quantify model robustness over time). 

#### 003_dynpars_utils.R

1. `get_dyn_pars` extracts dynamics parameters and plots model fits for each plot and spectral index. This function by default fits a Gompertz model with asymptotes constrained to 0 and 10. This is a wrapper for `nls_multstart::nls_multstart`.
2. `get_dyn_pars_lin` extracts dnyamics parameters for each plot and spectral index using linear interpolation. This function should be used if the temporal resolution of spectral measurements does not allow to fit parametric models. 
3. Various data wrangling helper functions called by 1. and 2. 

### Analysis

1. `01_prep_spcdat.R` Prepare spectral datasets for subquent analysis: (1) spectral index datasets, (2) spectral datasets for calibration and validation of (S)PLS-DA models. 
2. `02_pls_classification.R` Fit (S)PLS-DA models, extract variable importance, assess "internal", "general" and "dynamic" model performance. 
3. `03_extract_dynpars.R` Fit Gompertz model (optionally different models) and linear interpolation to spectral indices, extract model and dynamics parameters. 
4. `04_subset_selection.R` Select subsets of sensitive and insensitive spectral indices.
5. 

#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Fit gompertz model (optional: additional models) to SI data
# Extract model parameters from fits
# Extract dynamics parameters from fits

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#====================================================================================== -

.libPaths("")

library(tidyverse)
library(nlstools)
library(nls.multstart)
library(corrplot)

path_to_data <- ""
path_to_utils <- ""

source(paste0(path_to_utils, "003_dynpars_utils.R"))

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT/SI"))){
  dir.create(paste0(path_to_data, "OUT/SI"))
}
sinkdir <- paste0(path_to_data, "OUT/SI/") #Specify output directory

#====================================================================================== -

# FPWW023 =================== ----

#Load design, needed to create the final output df
inf <- readRDS(paste0(path_to_data, "prep_data/SIdat_FPWW023.rds")) %>% 
  dplyr::select(1:9) %>% unique() %>% mutate(Exp = "FPWW023")

#Load SI dataset 
DATA <- readRDS(paste0(path_to_data, "prep_data/SIdat_FPWW023.rds")) %>% 
  gather(SI, value, SI_760_730:SI_WI_NDVI_r)

#====================================================================================== -

#Fit Models ----
## Fit parametric models to SI time courses
#> all SI ----
## This is first done for all SI that seem amenable to this analysis
## This is evaluated graphically without model fits

ind_use <- c("SI_760_730", "SI_780_550", "SI_780_700",
             "SI_780_740", "SI_970_900_r","SI_ARVI", "SI_CARG",
             "SI_CARRE", "SI_CHLG","SI_CHLRE", 
             "SI_CIG","SI_CIRE", "SI_CLSI_r", 
             "SI_CNHI_ALI","SI_DCNI_ASD", "SI_DSWI",
             "SI_EVI", "SI_FII_r", "SI_GBNDVI_r", 
             "SI_GCC_HA", "SI_GM", "SI_GNDVI_HI",
             "SI_HI", "SI_LCI2", "SI_LWI",
             "SI_MCARI1", "SI_MCARI2", "SI_mND705",
             "SI_MSAVI", "SI_MSR_rev", "SI_MTCI", 
             "SI_MTVI1", "SI_MTVI2", "SI_NDMI", "SI_NDNI",
             "SI_NDRE", "SI_NDREI", "SI_NDRI_r",
             "SI_NDVI_nb_ASD", "SI_NDWI1", "SI_NDWI1650",
             "SI_NDWI2130", "SI_NGRDI", "SI_NHI_ALI", 
             "SI_NPCI_r", "SI_OSAVI", "SI_PBI",
             "SI_PRInorm_r", "SI_PSND1", "SI_PSND2",
             "SI_PSND3", "SI_PSND4", "SI_PSRI_r",
             "SI_PSSR1", "SI_PSSR2", "SI_PSSR3", 
             "SI_REIP", "SI_REM", "SI_RGR_r",
             "SI_RGR2_r", "SI_RVI1", "SI_SARVI", 
             "SI_SAVI", "SI_SIPI_r","SI_SIWSI",
             "SI_SLAIDI1", "SI_SR", "SI_SRWI", 
             "SI_TVI", 
             "SI_VI700", "SI_VIgreen", "SI_VIopt", 
             "SI_VOG1", "SI_VOG2_r", "SI_VOG3_r", 
             "SI_WDRVI", "SI_WI", "SI_YCAR")

#SI that must be fitted against thermal time excluding at least the last 3 measurement dates
dropfinal <- c("SI_DCNI_ASD", "SI_GBNDVI_r",
               "SI_HI", "SI_NDRI_r", "SI_NDSVI",
               "SI_NGRDI", "SI_REIP", "SI_RGR_r",
               "SI_RGR2_r", "SI_VIgreen")

Data <- DATA %>% filter(SI %in% ind_use)

#extract parameters for full phase
data <- Data %>% 
  filter(!SI %in% dropfinal)
pars <- split(data, data$SI) %>% lapply(get_dyn_pars, create_plot = TRUE)

#save pars data
if(!dir.exists(paste0(path_to_data, "OUT/SI/pars"))){
  dir.create(paste0(path_to_data, "OUT/SI/pars"))
}
saveRDS(pars, paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_drop0.rds"))

#extract parameters for restricted phase
dropdates <- as.Date(c("2018-07-09", "2018-07-12"))
data <- Data %>% 
  filter(!meas_date %in% dropdates) %>% droplevels() %>% 
  filter(SI %in% dropfinal)
pars <- split(data, data$SI) %>% lapply(get_dyn_pars, create_plot = TRUE)
saveRDS(pars, paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_drop2.rds"))

#====================================================================================== -

#> reduced SI ----
#SI set after exclusion of insatisfactorily modelled SI
ind_use <- c("SI_760_730", "SI_780_550", "SI_780_700",
             "SI_780_740", "SI_970_900_r","SI_ARVI",
             "SI_CHLRE", "SI_CIRE", "SI_CLSI_r", 
             "SI_DCNI_ASD", "SI_DSWI", "SI_EVI", 
             "SI_FII_r", "SI_GCC_HA", "SI_GNDVI_HI",
             "SI_HI", "SI_LCI2",
             "SI_MCARI1", "SI_MCARI2", "SI_mND705",
             "SI_MSAVI", "SI_MSR_rev", "SI_MTCI", 
             "SI_MTVI1", "SI_MTVI2", "SI_NDMI", "SI_NDNI",
             "SI_NDRE", "SI_NDREI", "SI_NDRI_r",
             "SI_NDVI_nb_ASD", "SI_NDWI1", "SI_NDWI1650",
             "SI_NDWI2130", "SI_NGRDI", "SI_NHI_ALI", 
             "SI_NPCI_r", "SI_OSAVI", "SI_PBI",
             "SI_PRInorm_r", "SI_PSND1", "SI_PSND2",
             "SI_PSND3", "SI_PSND4", "SI_PSRI_r",
             "SI_RGR2_r", "SI_SARVI", 
             "SI_SAVI", "SI_SIPI_r","SI_SIWSI",
             "SI_SLAIDI1", "SI_SR", "SI_SRWI", 
             "SI_VI700", "SI_VIgreen", "SI_VOG1", 
             "SI_WDRVI", "SI_WI", "SI_YCAR")

#SI that must be fitted against thermal time excluding at least the last 2 measurement dates
dropfinal <- c("SI_DCNI_ASD", "SI_NPCI_r",
               "SI_HI", "SI_NDRI_r", "SI_NDSVI",
               "SI_NGRDI", "SI_RGR2_r", "SI_VIgreen")

Data <- DATA %>% filter(SI %in% ind_use)

#for SI that can exploit the full phase
data <- Data %>% 
  filter(!SI %in% dropfinal)
pars <- split(data, data$SI) %>% lapply(get_dyn_pars, create_plot = TRUE)

#save pars data
if(!dir.exists(paste0(path_to_data, "OUT/SI/pars"))){
  dir.create(paste0(path_to_data, "OUT/SI/pars"))
}
saveRDS(pars, paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_drop0.rds"))

#for SI requiring exclusion of last 2 dates
dropdates <- as.Date(c("2018-07-09", "2018-07-12"))
data <- Data %>% 
  filter(!meas_date %in% dropdates) %>% droplevels() %>% 
  filter(SI %in% dropfinal)
pars <- split(data, data$SI) %>% lapply(get_dyn_pars, create_plot = TRUE)
saveRDS(pars, paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_drop2.rds"))

#====================================================================================== -

# Combine all data

#load data
data <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_drop0.rds"))
data2 <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_drop2.rds"))
data[names(data2)] <- data2

#drop NDRI and CLSI, dropped in second subset selection step
data <- data[c(-55, -9)]
saveRDS(data, paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_all.rds"))

########################################################################################################### -
########################################################################################################### -
########################################################################################################### -

#FPWW022 =================== ----
## The dynamics parameters are also extracted for the main experiment, 
## which will be used as indipendent validation 

#Load design, needed to create the final output df
inf <- readRDS(paste0(path_to_data, "prep_data/SIdat_FPWW022.rds")) %>% 
  dplyr::select(1:9) %>% 
  unique()

#Load SI dataset 
DATA <- readRDS(paste0(path_to_data, "prep_data/SIdat_FPWW022.rds")) %>% 
  gather(SI, value, SI_760_730:SI_WI_NDVI_r)

#Fit Models ----

#same SI as in FPWW023
ind_use_drop0 <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_drop0.rds")) %>% names()
ind_use_drop2 <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_drop2.rds")) %>% names()
ind_use <- c(ind_use_drop0, ind_use_drop2)

Data <- DATA %>% filter(SI %in% ind_use) %>% filter(!SI %in% c("SI_CLSI_r", "SI_NDRI_r"))

#for SI that can exploit the full phase
data <- Data %>% 
  filter(!SI %in% dropfinal)
pars <- split(data, data$SI) %>% lapply(get_dyn_pars, create_plot = FALSE)
saveRDS(pars, paste0(path_to_data, "OUT/SI/pars/pars_FPWW022_drop0.rds"))

#for SI requiring exclusion of last 2 dates
dropdates <- as.Date(c("2018-07-09", "2018-07-12"))
data <- Data %>% 
  filter(!meas_date %in% dropdates) %>% droplevels() %>% 
  filter(SI %in% dropfinal)
pars <- split(data, data$SI) %>% lapply(get_dyn_pars, create_plot = FALSE)
saveRDS(pars, paste0(path_to_data, "OUT/SI/pars/pars_FPWW022_drop2.rds"))

#====================================================================================== -

# Combine all data

#load data
data <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW022_drop0.rds"))
data2 <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW022_drop2.rds"))
data[names(data2)] <- data2

#drop NDRI and CLSI, dropped in second subset selection step
data <- rlist::list.remove(data, c("SI_NDRI_r", "SI_CLSI_r"))
saveRDS(data, paste0(path_to_data, "OUT/SI/pars/pars_FPWW022_all.rds"))

########################################################################################################### -
########################################################################################################### -
########################################################################################################### -

#FPWW012 =================== ----
## separate experiment carried out in 2016
## here, parameters are extracted using linear interpolation

#Load SI dataset 
data <- data.table::fread(paste0(path_to_data, "raw_data/SIdat_FPWW012.csv")) %>% as_tibble() %>% 
  dplyr::select(-c(grading_date:ref_date)) %>% 
  gather(SI, value, SI_760_730:SI_WI_NDVI_r) %>% group_by(Plot_ID)

#Load design, needed to create the final output df
inf <- data.table::fread(paste0(path_to_data, "raw_data/SIdat_FPWW012.csv")) %>% as_tibble() %>% 
  dplyr::select(Exp,Plot_ID, Lot,RangeLot, RowLot, Gen_Name, heading_date, heading_DAS, heading_GDDAS) %>%
  filter(complete.cases(.)) %>% 
  unique()

#same SI as in FPWW023
ind_use_drop0 <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_drop0.rds")) %>% names()
ind_use_drop2 <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_drop2.rds")) %>% names()
ind_use <- c(ind_use_drop0, ind_use_drop2)

#get the dynamics parameters for each SVI and each experimental plot
pars <- data %>% 
  filter(SI %in% ind_use) %>% 
  get_dyn_pars_lin(., create_plot = FALSE)

#drop NDRI and CLSI, dropped in second subset selection step
data <- rlist::list.remove(pars, c("SI_NDRI_r", "SI_CLSI_r"))
saveRDS(data, paste0(path_to_data, "OUT/SI/pars/pars_FPWW012_all.rds"))

#====================================================================================== -

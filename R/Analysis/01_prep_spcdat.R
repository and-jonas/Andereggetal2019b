#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Prepare spectral datasets for subsequent analyses
# 1) Spectral index dataset for FPWW023 and FPWW022
# 2) Spectral datasets for calibration and validation of a PLS-DA model

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

library(data.table)
library(tidyverse)
library(prospectr)
library(stringr)

path_to_utils <- "" #Specify path to utils folder
path_to_data <- "" #Specify path to data folder

#create sink directory
if(!dir.exists(paste0(path_to_data, "prep_data"))){
  dir.create(paste0(path_to_data, "prep_data"))
}
sinkdir <- paste0(path_to_data, "prep_data/") #Specify output directory

source(paste0(path_to_utils, "001_spectra_utils.R"))

#====================================================================================== -

#load designs, heading dates and temperature data
design_FPWW022 <- read.csv(paste0(path_to_data, "helper_files/design_FPWW022.csv")) %>% as_tibble()
design_FPWW023 <- read.csv(paste0(path_to_data, "helper_files/design_FPWW023.csv")) %>% as_tibble()
heading <- read.csv(paste0(path_to_data, "raw_data/heading.csv"), sep = ";") %>% as_tibble() %>% 
  mutate(heading_date = as.Date(heading_date, "%d.%m.%Y"))
temp <- read.csv(paste0(path_to_data, "raw_data/temp.csv")) %>% as_tibble()

#load spectral data
spc_FPWW023 <- data.table::fread(paste0(path_to_data, "raw_data/spectra_FPWW023_full.csv")) %>% 
  as_tibble() %>% dplyr::select(-spc_ID) %>% 
  mutate_at(vars(Plot_ID), as.factor) %>% 
  mutate(meas_date = as.Date(meas_date),
         replicate = as.numeric(replicate))
spc_FPWW022 <- data.table::fread(paste0(path_to_data, "raw_data/spectra_FPWW022.csv")) %>% 
  as_tibble() %>% dplyr::select(-spc_ID) %>% 
  mutate_at(vars(Plot_ID), as.factor) %>% 
  mutate(meas_date = as.Date(meas_date),
         replicate = as.numeric(replicate))

#====================================================================================== -

#> SPECTRAL INDEX DATASET ----
#>> FPWW023 (inoculated and corresponding control plots) ----

#Calculate spectral indices
SI <- spc_FPWW023 %>%
  f_spc_smooth(p = 3, w = 11, m = 0) %>%
  f_spc_avg() %>%
  f_calc_si() %>%
  f_scale_si()

#add design, heading, and measurement time point
fpww23 <- left_join(design_FPWW023, SI, by = "Plot_ID") %>%
  left_join(heading, by = "Gen_Name") %>%
  #rearrange columns
  dplyr::select(1:6, heading_date:heading_GDDAS, meas_date, everything()) %>% 
  #add treatment information
  mutate(trt = ifelse(grepl("FPWW023", Plot_ID), "dis", "ctrl"))
  
  #Calcualte meas_gddah as meas_GDD - heading_GDD for each plot and each meas date
  meas_GDDAH <- NULL
  for (i in 1:nrow(fpww23)) {
    if (!is.na(fpww23$heading_date)[i]) {
      meas_GDDAH[i] <- temp[temp$day == paste(fpww23$meas_date)[i], 9] -
        temp[temp$day == paste(fpww23$heading_date)[i], 9]
    }
    else paste("NA")
    print(meas_GDDAH[i])
  }
  fpww23$meas_GDDAH <- unlist(meas_GDDAH)

fpww23 <- fpww23 %>% 
  #normalize meas_GDDAH for differences in heading date
  group_by(Plot_ID) %>%
  mutate(meas_GDDAH_norm = meas_GDDAH - first(meas_GDDAH)) %>%
  ungroup() %>% 
  #rearrange columns
  dplyr::select(1:10, meas_GDDAH, meas_GDDAH_norm, trt, everything())
saveRDS(fpww23, paste0(sinkdir, "SIdat_FPWW023.rds"))

#====================================================================================== -

#>> FPWW022 (main experiment) ----

#Calculate spectral indices
SI <- spc_FPWW022 %>%
  f_spc_smooth(p = 3, w = 11, m = 0) %>%
  f_spc_avg() %>%
  f_calc_si() %>%
  f_scale_si()

#add design, heading, and measurement time point
fpww22 <- right_join(design_FPWW022, SI, by = "Plot_ID") %>%
  left_join(heading, by = "Gen_Name") %>%
  dplyr::select(1:6, heading_date:heading_GDDAS, meas_date, everything()) %>% 
  mutate(trt = ifelse(grepl("FPWW023", Plot_ID), "dis", "ctrl"))

  #Calcualte meas_gddah as meas_GDD - heading_GDD for each plot and each meas date
  meas_GDDAH <- NULL
  for (i in 1:nrow(fpww22)) {
    if (!is.na(fpww22$heading_date)[i]) {
      meas_GDDAH[i] <- temp[temp$day == paste(fpww22$meas_date)[i], 9] -
        temp[temp$day == paste(fpww22$heading_date)[i], 9]
    }
    else paste("NA")
    print(i)
  }
  fpww22$meas_GDDAH <- unlist(meas_GDDAH)

fpww22 <- fpww22 %>% 
  #normalize meas_GDDAH for differences in heading date
  group_by(Plot_ID) %>%
  mutate(meas_GDDAH_norm = meas_GDDAH - first(meas_GDDAH)) %>% 
  ungroup() %>% 
  #rearrange columns
  dplyr::select(1:10, meas_GDDAH, meas_GDDAH_norm, trt, everything())
saveRDS(fpww22, paste0(sinkdir, "SIdat_FPWW022.rds"))

#====================================================================================== -

#> FULL-SPECTRUM DATASETS ----
#>> FPWW023 ----

#Pre-process spectral data
#bin spectral data with bin size = 3
#add info
spc_pls <- spc_FPWW023 %>%
  f_spc_smooth(p = 3, w = 11, m = 0) %>%
  f_spc_avg() %>%
  f_spc_bin(bin_size = 3) %>% 
  f_spc_trim() %>% 
  #add plot and measurement information
  left_join(fpww23[1:13], by = c("Plot_ID", "meas_date")) %>% 
  #rearrange columns
  dplyr::select(Plot_ID, meas_date, Exp:trt, everything())
saveRDS(spc_pls, paste0(sinkdir, "spectra_plsda_bin3.rds"))

#Pre-process spectral data
#bin spectral data with bin size = 6
#add info
spc_pls_fpww23 <- spc_FPWW023 %>%
  f_spc_smooth(p = 3, w = 11, m = 0) %>%
  f_spc_avg() %>%
  f_spc_bin(bin_size = 6) %>% 
  f_spc_trim() %>% 
  left_join(fpww23[1:13], by = c("Plot_ID", "meas_date")) %>% 
  dplyr::select(Plot_ID, meas_date, Exp:trt, everything())
saveRDS(spc_pls_fpww23, paste0(sinkdir, "spectra_plsda_bin6.rds"))

#>> FPWW022 ----
spc_pls_fpww22 <- spc_FPWW022 %>%
  f_spc_smooth(p = 3, w = 11, m = 0) %>%
  f_spc_avg() %>%
  f_spc_bin(bin_size = 6) %>%
  f_spc_trim()
saveRDS(spc_pls_fpww22, paste0(sinkdir, "spectra_plsda_FPWW022_bin6.rds"))

#====================================================================================== -

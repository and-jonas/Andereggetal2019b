#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Select subset of sensitive and insensitive SI

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
library(corrplot)

path_to_data <- ""

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT/SI"))){
  dir.create(paste0(path_to_data, "OUT/SI"))
}
sinkdir <- paste0(path_to_data, "OUT/SI/") #Specify output directory

#====================================================================================== -

data <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_all.rds"))

# Select pot. normalizers ----
#select SI with minimal shift in "checkpoints"
keypoints_norm <- do.call("rbind", data) %>% tibble::as_tibble() %>% 
  group_by(trt, SI, par) %>% 
  summarise(mean = mean(val)) %>% 
  tidyr::spread(trt, mean) %>% 
  dplyr::filter(par %in% c("M", "t85")) %>% 
  mutate(diff = ctrl-dis) %>% 
  arrange(par, diff) %>% 
  group_by(par) %>% 
  #top 8 for each parameter
  top_n(n = 8, wt = desc(abs(diff))) %>% ungroup() %>% 
  distinct(SI) %>% pull(SI)

#select SI with minimal change in duration
dur_norm <- do.call("rbind", data) %>% tibble::as_tibble() %>% 
  group_by(trt, SI, par) %>% 
  summarise(mean = mean(val, na.rm = TRUE)) %>% 
  tidyr::spread(trt, mean) %>% 
  dplyr::filter(par %in% c("dur1")) %>% 
  mutate(rat = dis/ctrl) %>% 
  arrange(par, desc(rat)) %>% 
  group_by(par) %>% 
  #top 8 for each parameter
  top_n(n = 8, wt = abs(rat)) %>% ungroup() %>% 
  distinct(SI) %>% pull(SI)

#select SI with minimal change in rate
rat_norm <- do.call("rbind", data) %>% tibble::as_tibble() %>% 
  group_by(trt, SI, par) %>% 
  summarise(mean = mean(val, na.rm = TRUE)) %>% 
  tidyr::spread(trt, mean) %>% 
  dplyr::filter(par == "b") %>% 
  mutate(rat = dis/ctrl) %>% 
  arrange(par, desc(rat)) %>% 
  group_by(par) %>% 
  #top 8 for each parameter
  top_n(n = 8, wt = abs(rat)) %>% ungroup() %>% 
  distinct(SI) %>% pull(SI)

#====================================================================================== -
#====================================================================================== -
  
# Select keypoints norm subset ----
tidy <- data %>% 
  do.call("rbind", .) %>% 
  filter(par %in% c("t85", "M")) %>% 
  mutate(par = paste(SI, par, sep = "-")) %>% 
  group_by(SI) %>% nest() %>% 
  mutate(data = purrr::map(.x = data, .f = spread, par, val) %>% 
           purrr::map(., .f = dplyr::select, starts_with("SI_")))

norm <- tidy %>% 
  filter(SI %in% keypoints_norm) %>% 
  spread(SI, data) %>% 
  unnest()

normlist <- list(norm[,grepl("-M", names(norm))], 
                 norm[,grepl("-t85", names(norm))])

#drop redundant SI iterativelly
sel_normsi <- lapply(normlist, function(x){
  x <- dplyr::select(x, -matches("780_550|PSND4|780_700|PBI|MTCI|PSND3"))
})

# #Select subset based on pair-wise correlations
pdf(paste(paste0(path_to_data, "OUT/SI/pars/corr_norm_kp.pdf"), sep = ""), 7, 5)
for(i in 1:length(sel_normsi)){
  d <- sel_normsi[[i]]
  M <- cor(d, use = "pairwise.complete.obs", method = "pearson")
  summary_stats <- list(mean(M[lower.tri(M)]), min(M[lower.tri(M)]), max(M[lower.tri(M)]))
  corrplot(M, method = "color", tl.cex = 0.5, tl.col = "black", addCoef.col = "black", number.cex = 0.5)
}
dev.off()

#====================================================================================== -

# Select keypoints diff subset ----

diff <- tidy %>% 
  filter(!SI %in% keypoints_norm) %>% 
  spread(SI, data) %>% 
  unnest() 

difflist <- list(diff[,grepl("-M", names(diff))], 
                 diff[,grepl("-t85", names(diff))])

#drop redundant SI iteratively
sel_diffsi <- lapply(difflist, function(x){
  x <- dplyr::select(x, -matches("MCARI1|CIRE-|ARVI-|EVI|MTVI1|MTVI2|760_730|NDMI|NDRE-|NDREI???|LCI2|NDWI1650|NDWI1-|mnD705|GCC_HA|NPCI|PSND1|PSND2|SI_OSAVI-|RGR2|MSAVI|SIWSI|SRWI|SLAIDI|VIgreen|MSR_rev|WDRVI|970|NHI_ALI|SAVI"))
})

# #Select subset based on pair-wise correlations
pdf(paste0(path_to_data, "OUT/SI/pars/corr_diff_kp.pdf"), 7, 5)
for(i in 1:length(sel_diffsi)){
  d <- sel_diffsi[[i]]
  M <- cor(d, use = "pairwise.complete.obs", method = "pearson")
  summary_stats <- list(mean(M[lower.tri(M)]), min(M[lower.tri(M)]), max(M[lower.tri(M)]))
  corrplot(M, method = "color", tl.cex = 0.5, tl.col = "black", addCoef.col = "black", number.cex = 0.5)
}
dev.off()

#====================================================================================== -

#save selected SI
if(!dir.exists(paste0(path_to_data, "OUT/SI/preds"))){
  dir.create(paste0(path_to_data, "OUT/SI/preds"))
}
lapply(str_split(names(sel_normsi[[1]]), "-"), "[[", 1) %>% unlist() %>% 
  saveRDS(., paste0(path_to_data, "OUT/SI/preds/normsi_kp.rds"))
lapply(str_split(names(sel_diffsi[[1]]), "-"), "[[", 1) %>% unlist() %>% 
  saveRDS(., paste0(path_to_data, "OUT/SI/preds/diffsi_kp.rds"))

#====================================================================================== -

# Select ratepars norm subset ----

tidy <- data %>% 
  do.call("rbind", .) %>% 
  filter(par %in% c("b", "dur1")) %>% 
  mutate(par = paste(SI, par, sep = "-")) %>% 
  group_by(SI) %>% nest() %>% 
  mutate(data = purrr::map(.x = data, .f = spread, par, val) %>% 
           purrr::map(., .f = dplyr::select, starts_with("SI_")))

norm <- tidy %>% 
  filter(SI %in% rat_norm) %>% 
  spread(SI, data) %>% 
  unnest()

normlist <- list(norm[,grepl("-b", names(norm))], 
                 norm[,grepl("-dur1", names(norm))])

#drop redundant SI iterativelly
sel_normsi <- lapply(normlist, function(x){
  x <- dplyr::select(x, -matches("NDNI|PBI|780_700|MSR_rev"))
})

# #Select subset based on pair-wise correlations
pdf(paste0(path_to_data, "OUT/SI/pars/corr_norm_rat.pdf"), 7, 5)
for(i in 1:length(sel_normsi)){
  d <- sel_normsi[[i]]
  M <- cor(d, use = "pairwise.complete.obs", method = "pearson")
  summary_stats <- list(mean(M[lower.tri(M)]), min(M[lower.tri(M)]), max(M[lower.tri(M)]))
  corrplot(M, method = "color", tl.cex = 0.65, tl.col = "black", addCoef.col = "black")
}
dev.off()

#====================================================================================== -

# Select ratepars diff subset ----

diff <- tidy %>% 
  filter(!SI %in% keypoints_norm) %>% 
  spread(SI, data) %>% 
  unnest() 

difflist <- list(diff[,grepl("-b", names(diff))], 
                 diff[,grepl("-dur1", names(diff))])

#drop the same as for keypoints;
#drop additional ones if necessary
diffsi <- readRDS(paste0(path_to_data, "OUT/SI/preds/diffsi_kp.rds"))
sel_diffsi <- lapply(difflist, function(x){
  x <- dplyr::select(x, matches(paste(diffsi, collapse = "|")), -matches("SR-|NDNI"))
})

#Select subset based on pair-wise correlations
pdf(paste0(path_to_data, "OUT/SI/pars/corr_diff_rat.pdf"), 7, 5)
for(i in 1:length(sel_diffsi)){
  d <- sel_diffsi[[i]]
  M <- cor(d, use = "pairwise.complete.obs", method = "pearson")
  summary_stats <- list(mean(M[lower.tri(M)]), min(M[lower.tri(M)]), max(M[lower.tri(M)]))
  corrplot(M, method = "color", tl.cex = 0.5, tl.col = "black", addCoef.col = "black", number.cex = 0.5)
}
dev.off()

#====================================================================================== -

#save selected SI
lapply(str_split(names(sel_normsi[[1]]), "-"), "[[", 1) %>% unlist() %>% 
  saveRDS(., paste0(path_to_data, "OUT/SI/preds/normsi_rat.rds"))
lapply(str_split(names(sel_diffsi[[1]]), "-"), "[[", 1) %>% unlist() %>% 
  saveRDS(., paste0(path_to_data, "OUT/SI/preds/diffsi_rat.rds"))

#====================================================================================== -
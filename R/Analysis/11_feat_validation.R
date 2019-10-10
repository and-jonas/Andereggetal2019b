#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Validation of spectral-temporal features on a 2016 experiment. 

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

source(paste0(path_to_utils, "003_dynpars_utils.R"))

#load spectral index dataset
data <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW012_all.rds")) %>% split(., .$SI)

#load STB dat
placl <- read.csv(paste0(path_to_data, "raw_data/stb_FPWW012.csv")) %>% 
  spread(timestamp, mean_placl) %>% 
  rename(mean_placl_20160520 = 2,
         mean_placl_20160704 = 3)

#identify plots affected by lodging and ggt
lodg_plots <- read.csv(paste0(path_to_data, "raw_data/lodging.csv"), sep = ";") %>%
  filter(Lodging %in% c("l", "ll")) %>% pull(Plot_ID) %>% as.character()
ggt_plots <- read.csv(paste0(path_to_data, "raw_data/GGT_Rating.csv"), skip = 3, sep = ";") %>% 
  arrange(Plot_ID) %>% as_tibble() %>% rename("ggt" = 3) %>% 
  filter(ggt > 1) %>% pull(Plot_ID) %>% as.character()
rm_plots <- c(lodg_plots, ggt_plots) %>% unique()

#====================================================================================== -

#define selected spectral indices for KEYPOINTS
normsi <- readRDS(paste0(path_to_data, "OUT/SI/preds/normsi_kp.rds"))
diffsi <- readRDS(paste0(path_to_data, "OUT/SI/preds/diffsi_kp.rds"))
#select data for normalizers and differentiators
dat_norm <- data[normsi] %>%  
  lapply(., select_pars, which = "norm", pars = c("t85", "t50"))
dat_diff <- data[diffsi] %>% 
  lapply(., select_pars, which = "diff", pars = c("t85", "t50"))
#extract DELTAS
deltas <- list()
for(i in normsi){
  deltas[[i]] <- dat_diff %>% 
    #for each diff SI
    lapply(function(x) x %>%
             full_join(dat_norm[[i]]) %>% 
             mutate(t50_delta = t50_diffind - t50_normind,
                    t85_delta = t85_diffind - t85_normind) %>% 
             dplyr::select(-contains("diffind"), -contains("normind")) %>% 
             tidyr::gather(par, delta, t50_delta:t85_delta))
  deltas[[i]] <- do.call("rbind", deltas[[i]])
}
deltas <- do.call("rbind", deltas)

del_mat <- deltas %>% 
  mutate(pred = paste(SI_diff, SI_norm, par, sep = "-")) %>% 
  dplyr::select(-SI_diff, -SI_norm, -par) %>% 
  tidyr::spread(pred, delta)

# add response variable

preds <- full_join(del_mat, placl, by = "Plot_ID") %>% 
  dplyr::select(Plot_ID, mean_placl_20160520, mean_placl_20160704, contains("SI_")) %>% 
  filter(!Plot_ID %in% rm_plots) %>% droplevels()

#====================================================================================== -

#Correlation analyses ----
feats <- preds %>% dplyr::select(contains("SI_")) %>% names()
cor <- vector()
for(i in feats){
  cor[i] <- abs(cor(preds$mean_placl_20160704, preds[i], use = "pairwise.complete")) %>% 
    as.data.frame(.) %>% round(., 2)
}
cors <- unlist(cor) %>% as.data.frame() %>% 
  tibble::rownames_to_column() %>% as_tibble() %>% 
  rename(., Pearson_r = .) %>% 
  arrange(desc(Pearson_r))
ggplot(preds) +
  geom_point(aes(x =abs(`SI_MCARI2-SI_SIPI_r-t50_delta`), y = mean_placl_20160704), shape = 16, size = 3.25, alpha = 0.4) +
  xlab("|MCARI2-SIPI-t50_delta|") + ylab("PLACL") +
  geom_smooth(aes(x = abs(`SI_MCARI2-SI_SIPI_r-t50_delta`), y = mean_placl_20160704), method = "lm", color = "blue", size = 1.25) +
  annotate(x = 50, y = 75, "text", label = "r = 0.51") +
  theme_bw() +
  theme(panel.grid = element_blank())

#====================================================================================== -
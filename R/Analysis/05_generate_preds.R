#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Generate spectral-temporal features (deltas and ratios)
# Assemble predictor matrix

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

path_to_data <- ""
path_to_utils <- ""

source(paste0(path_to_utils, "003_dynpars_utils.R"))

FPWW023 <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW023_all.rds"))
FPWW022 <- readRDS(paste0(path_to_data, "OUT/SI/pars/pars_FPWW022_all.rds"))
data <- list("FPWW023" = FPWW023, 
             "FPWW022" = FPWW022)

create_plots <- list("FPWW023" = TRUE, 
                     "FPWW022" = FALSE)

#====================================================================================== -

#for each of the Experiments FPWW022 and FPWW023
#extract spectral-temporal features
for(k in names(data)){
  #Extract deltas and ratios ----
  #> Deltas ----
  #define selected spectral indices for KEYPOINTS
  normsi <- readRDS(paste0(path_to_data, "OUT/SI/preds/normsi_kp.rds"))
  diffsi <- readRDS(paste0(path_to_data, "OUT/SI/preds/diffsi_kp.rds"))
  #select data for normalizers and differentiators
  dat_norm <- data[[k]][normsi] %>%  
    lapply(., select_pars, which = "norm", pars = c("t85", "M"))
  dat_diff <- data[[k]][diffsi] %>% 
    lapply(., select_pars, which = "diff", pars = c("t85", "M"))
  #extract DELTAS
  deltas <- list()
  for(i in normsi){
    deltas[[i]] <- dat_diff %>% 
      #for each diff SI
      lapply(function(x) x %>%
               full_join(dat_norm[[i]]) %>% 
               mutate(M_delta = M_diffind - M_normind,
                      t85_delta = t85_diffind - t85_normind) %>% 
               dplyr::select(-contains("diffind"), -contains("normind")) %>% 
               tidyr::gather(par, delta, M_delta:t85_delta))
    deltas[[i]] <- do.call("rbind", deltas[[i]])
  }
  deltas <- do.call("rbind", deltas)
  #====================================================================================== -
  #> Ratios ----
  #define selected spectral indices for CHANGE PARAMETERS
  normsi <- readRDS(paste0(path_to_data, "OUT/SI/preds/normsi_rat.rds"))
  diffsi <- readRDS(paste0(path_to_data, "OUT/SI/preds/diffsi_rat.rds"))
  #select data for normalizers and differentiators
  dat_norm <- data[[k]][normsi] %>%  
    lapply(., select_pars, which = "norm", pars = c("b", "dur1"))
  dat_diff <- data[[k]][diffsi] %>% 
    lapply(., select_pars, which = "diff", pars = c("b", "dur1"))
  #extract RATIOS
  ratios <- list()
  for(i in normsi){
    ratios[[i]] <- dat_diff %>% 
      #for each diff SI
      lapply(function(x) x %>%
               full_join(dat_norm[[i]]) %>% 
               mutate(b_ratio = b_diffind/b_normind,
                      dur1_ratio = dur1_diffind/dur1_normind) %>% 
               dplyr::select(-contains("diffind"), -contains("normind")) %>% 
               tidyr::gather(par, ratio, b_ratio:dur1_ratio))
    ratios[[i]] <- do.call("rbind", ratios[[i]])
  }
  ratios <- do.call("rbind", ratios)
  #====================================================================================== -
  # Create boxplots of deltas and ratios ----
  plot <- create_plots[[k]]
  if(plot){
    #> Deltas ----
    pdf(paste0(path_to_data, "OUT/SI/preds/deltas.pdf"), 7, 5)
    #for each normalization SI
    for(i in unique(deltas$SI_norm)){
      #select data
      df <- deltas %>% dplyr::filter(SI_norm == i)
      #for each differentiation SI
      for(j in unique(df$SI_diff)){
        ddf <- df %>% filter(!Plot_ID == "FPWW0220320" & SI_diff == j)
        plot <- ggplot(ddf) +
          geom_boxplot(aes(x = trt, y = delta)) +
          geom_jitter(aes(x = trt, y = delta, col = Lot),
                      position=position_jitter(width=0.1,height=0),
                      alpha=0.5,
                      size=1.5,
                      show.legend = TRUE) +
          theme_bw() +
          facet_wrap(~par, scales = "free") +
          ggtitle(paste("diff_ind = ", j, "\nnorm_ind = ", i, sep = ""))
        print(plot)
      }
    }
    dev.off()
    #> Ratios ----
    pdf(paste0(path_to_data, "OUT/SI/preds/ratios.pdf"), 7, 5)
    #for each normalization SI
    for(i in unique(ratios$SI_norm)){
      df <- ratios %>% dplyr::filter(SI_norm == i)
      for(j in unique(df$SI_diff)){
        ddf <- df %>% filter(!Plot_ID == "FPWW0220320" & SI_diff == j)
        plot <- ggplot(ddf) +
          geom_boxplot(aes(x = trt, y = ratio)) +
          geom_jitter(aes(x = trt, y = ratio, col = Lot),
                      position=position_jitter(width=0.1,height=0),
                      alpha=0.5,
                      size=1.5,
                      show.legend = TRUE) +
          theme_bw() +
          facet_wrap(~par, scales = "free") +
          ggtitle(paste("diff_ind = ", j, "\nnorm_ind = ", i, sep = ""))
        print(plot)
      }
    }
    dev.off()
  }
  #====================================================================================== -
  # Assemble predictor matrix ----
  rat_mat <- ratios %>% 
    mutate(pred = paste(SI_diff, SI_norm, par, sep = "-")) %>% 
    dplyr::select(-SI_diff, -SI_norm, -par) %>% 
    tidyr::spread(pred, ratio) 
  del_mat <- deltas %>% 
    mutate(pred = paste(SI_diff, SI_norm, par, sep = "-")) %>% 
    dplyr::select(-SI_diff, -SI_norm, -par) %>% 
    tidyr::spread(pred, delta)
  pred_mat <- full_join(del_mat, rat_mat, by = c("Plot_ID", "Exp", "Lot", "RangeLot", "RowLot",
                                                 "Gen_Name", "trt", "heading_date", "heading_DAS", 
                                                 "heading_GDDAS"))
  saveRDS(pred_mat, paste0(path_to_data, "OUT/SI/preds/predmat", k, ".rds"))
}

#====================================================================================== -

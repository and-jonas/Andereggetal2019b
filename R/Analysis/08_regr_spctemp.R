#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Train regression models for disease severity prediction
# Test performance on FPWW022

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
library(caret)
library(ranger)
library(desplot)
library(Cubist)
library(elasticnet)
library(pls)

path_to_data <- ""
path_to_utils <- ""

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT/Models_spctemp/regr"))){
  dir.create(paste0(path_to_data, "OUT/Models_spctemp/regr"))
}
sinkdir <- paste0(path_to_data, "OUT/Models_spctemp/regr/") #Specify output directory

source(paste0(path_to_utils, "004_plot_predobs.R"))

#====================================================================================== -

#Prepare data ----

#get predictor data
spc0 <- readRDS(paste0(path_to_data, "OUT/SI/preds/predmatFPWW023.rds"))

#> (1) use all controls
spc <- spc0

#get disease severity data
dis_sev <- read.csv(paste0(path_to_data, "raw_data/stb_sev.csv")) %>% dplyr::select(Plot_ID, sev)
DATA <- left_join(spc, dis_sev)

#record plot sequence
design <- left_join(spc, dis_sev) %>% dplyr::select(-contains("SI_"), -trt, -sev, -contains("heading"))
saveRDS(design, paste0(path_to_data, "helper_files/design_feat_all.rds"))

#prepare data for modelling
data <- DATA %>% 
  select_if(function(x) !any(is.na(x))) %>% 
  mutate(trt = as.factor(trt)) %>% tidyr::drop_na() %>% 
  dplyr::select(sev, contains("SI_"))
#transform varnames to "legal" names
names(data) <- c("sev", paste("V", 1:(length(data)-1), sep = ""))
saveRDS(data, paste0(path_to_data, "OUT/SI/preds/predmat_regr_allctrl.rds"))

#====================================================================================== -

#> (2) use only inoculated plots
spc <- spc0 %>% dplyr::filter(trt == "dis")
#record plot sequence
design <- left_join(spc, dis_sev) %>% dplyr::select(-contains("SI_"), -trt, -sev, -contains("heading"))
saveRDS(design, paste0(path_to_data, "helper_files/design_feat_sub.rds"))

#prepare data for modelling
data <- DATA %>% 
  select_if(function(x) !any(is.na(x))) %>% 
  mutate(trt = as.factor(trt)) %>% tidyr::drop_na() %>% 
  dplyr::select(sev, contains("SI_"))
#transform varnames to "legal" names
names(data) <- c("sev", paste("V", 1:(length(data)-1), sep = ""))
saveRDS(data, paste0(path_to_data, "OUT/SI/preds/predmat_regr_noctrl.rds"))

#====================================================================================== -

data <- as.list(list.files(paste0(path_to_data, "OUT/SI/preds/"), pattern = "predmat_regr_", full.names = T)) %>% 
  lapply(., readRDS)
names(data) <- list("allctrl", "noctrl")

# FULL MODELS =============== ----
## for the full dataset and a dataset only comprising the inoculated plots
## evaluate different models and save result

for (i in names(data)){
  
  #get data
  dat <- data[[i]]
  #> PLSR ----
  indx <- createMultiFolds(dat$sev, k = 10, times = 10)
  ctrl <- caret::trainControl(method = "repeatedcv", 
                              index = indx,
                              savePredictions = TRUE,
                              verboseIter = TRUE,
                              selectionFunction = "oneSE")
  plsr <- caret::train(sev ~ .,
                       data = dat,
                       preProcess = c("center", "scale"),
                       method = "pls",
                       tuneLength = 15, 
                       trControl = ctrl,
                       importance = TRUE)
  saveRDS(plsr, paste0(sinkdir, "plsr_", i, ".rds"))
  #> Cubist ----
  #grid of values to test
  train_grid <- expand.grid(committees = c(1, 2, 3, 5, 10, 15, 20, 50),
                            neighbors = 0)
  cubist <- train(sev ~ .,
                  data = dat,
                  preProcess = c("center", "scale"),
                  method = "cubist",
                  tuneGrid = train_grid,
                  trControl = ctrl,
                  importance = TRUE)
  saveRDS(cubist, paste0(sinkdir, "cubist_", i, ".rds"))
  #> rf-ranger ----
  #grid of mtry values to test
  mtry <- c(1, 2, 5, 9, 14, 20, 30, 45, 70, 100, 200)
  min_nodes <- c(1, 2, 5, 10)
  #specify model tuning parameters
  tune_grid <- expand.grid(mtry = mtry,
                           splitrule = "variance", #default
                           min.node.size = min_nodes)
  rf_ranger <- caret::train(sev ~ .,
                            data = dat,
                            preProc = c("center", "scale"),
                            method = "ranger",
                            tuneGrid = tune_grid,
                            importance = "permutation",
                            num.trees = 2000,
                            num.threads = 10, 
                            trControl = ctrl)
  saveRDS(rf_ranger, paste0(sinkdir, "rf_", i, ".rds"))
  #> Ridge ----
  ridgeGrid <- data.frame(lambda = seq(0, .15, length = 15))
  ridge <- caret::train(sev ~ .,
                        data = dat,
                        # preProcess = NULL,
                        preProcess = c("center", "scale"),
                        method = "ridge",
                        tuneGrid = ridgeGrid,
                        trControl = ctrl,
                        importance = TRUE)
  saveRDS(ridge, paste0(sinkdir, "ridge_", i, ".rds"))
  
}

#====================================================================================== -

#Create predobs plots ----

setwd(sinkdir)

mod_names <- as.list(list.files(pattern = "noctrl"))
mod_names_adj <- as.list(list.files(pattern = "allctrl"))

mods <- lapply(mod_names, readRDS)
mods_adj <- lapply(mod_names_adj, readRDS)

plots_all <- readRDS(paste0(path_to_data, "helper_files/design_feat_all.rds"))
plots_sub <- readRDS(paste0(path_to_data, "helper_files/design_feat_sub.rds"))

perf_all <- lapply(mods_adj, assess_mod_perf, plots_all)
figs_all <- lapply(perf_all, plot_predobs, adjust = TRUE)

perf_sub <- lapply(mods, assess_mod_perf, plots_sub)
figs_sub <- lapply(perf_sub, plot_predobs, adjust = FALSE)

#====================================================================================== -

#Model generalization ----

#plots included in the training set;
#never include these, to garantee same testing set for all approaches
train_plots <- spc0 %>% dplyr::filter(grepl("FPWW022", Plot_ID)) %>% pull(Plot_ID) %>% unique()
preds <- readRDS(paste0(path_to-data, "OUT/SI/preds/predmatFPWW022.rds"))

#prepare FPWW023 data
newdata <- preds %>% 
  #drop plots used for training
  dplyr::filter(!Plot_ID %in% train_plots) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(contains("SI_"))
names(newdata) <- paste("V", 1:(length(newdata)), sep = "")

plot_seq <- preds  %>% 
  #drop plots used for training
  dplyr::filter(!Plot_ID %in% train_plots) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  pull(Plot_ID)

#get plot information for desplot
plot_inf <- preds %>% 
  #drop plots used for training
  mutate(trt = as.factor(trt)) %>% 
  dplyr::select(-contains("SI_")) %>%
  dplyr::filter(Plot_ID %in% plot_seq)

#get severity predictons
preds <- lapply(mods_adj, predict.train, newdata = newdata[ , names(newdata) != "sev"]) %>% 
  lapply(., as.numeric) %>% 
  lapply(., expss::na_if, expss::gt(0.5))

plot_distr <- list()
pred_data <- list()
for(i in 1:length(preds)){
  pred_data[[i]] <- preds[i] %>% as.data.frame(col.names = "V1") %>% 
    mutate(model = paste(mod_names_adj[i]) %>% gsub(".rds", "", .))
  plot_distr[[i]] <- ggplot(pred_data[[i]]) +
    geom_bar(aes(x = V1), stat = "bin", binwidth = 0.01) +
    facet_wrap(~model) +
    scale_x_continuous(limits = c(-0.5, 0.5)) +
    scale_y_continuous(limits = c(0, 50)) +
    geom_vline(xintercept = 0, col = "red") +
    geom_vline(xintercept = 0.05, col = "red", lty = 2) +
    xlab("Predicted STB severity") + ylab("Count") +
    theme_bw()
}

#====================================================================================== -

#Spatial Patterns ----

pred_des <- lapply(pred_data, function(x) cbind(x, plot_seq) %>% 
  full_join(., plot_inf, by = c("plot_seq" = "Plot_ID")) %>% 
  as_tibble() %>% 
  rename(pred = 1, 
         algorithm = 2) %>% 
  mutate(pred = ifelse(pred < 0, 0, pred)))

plot <- list()
for(i in 1:length(pred_des)){
  data <- pred_des[[i]]
  plot[[i]] <- desplot(pred ~ RangeLot + RowLot | Lot,
                       data = data, cex = 0.6, ticks = TRUE, show.key = TRUE,
                       midpoint = "midrange",
                       col.regions = RedGrayBlue,
                       main = paste(unique(pred_des[[i]]$algorithm)))
}

#====================================================================================== -
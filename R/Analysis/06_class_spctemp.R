#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Train classifiers for diseased/disease-free 
# Test performance of classifiers on FPWW022
# Perform recursive feature elimination with rf as base-learner

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
library(stringr)
library(caret)
library(ranger)
library(desplot)

path_to_data <- ""

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT/Models_spctemp/class"))){
  dir.create(paste0(path_to_data, "OUT/Models_spctemp/class"), recursive = TRUE)
}
sinkdir <- paste0(path_to_data, "OUT/Models_spctemp/class/") #Specify output directory

#====================================================================================== -

#Prepare data ----

##Prepare "training data"
data_FPWW023 <- readRDS(paste0(path_to_data, "OUT/SI/preds/predmatFPWW023.rds"))
#extract variable names
var_names <- data_FPWW023 %>% dplyr::select(contains("SI_")) %>% names() %>% data.frame() #these must be dropped ("not legal names")
traindat <- data_FPWW023 %>% 
  select_if(function(x) !any(is.na(x))) %>% 
  mutate(trt = as.factor(trt)) %>% tidyr::drop_na() %>% 
  dplyr::select(trt, contains("SI_"))
names(traindat) <- c("trt", paste("V", 1:(length(traindat)-1), sep = ""))
saveRDS(traindat, paste0(path_to_data, "OUT/SI/preds/predmat_class.rds"))

#Prepare "testing data"
data_FPWW022 <- readRDS(paste0(path_to_data, "OUT/SI/preds/predmatFPWW022.rds")) 
testdat <- data_FPWW022 %>% #prepare FPWW023 for validation
  #drop plots used for training
  dplyr::filter(!Plot_ID %in% data_FPWW023$Plot_ID) %>% 
  dplyr::select(trt, contains("SI_")) %>% 
  mutate(trt = as.factor(trt)) %>% drop_na()
names(testdat) <- c("trt", paste("V", 1:(length(testdat)-1), sep = "")) 
stopifnot(length(traindat) == length(testdat))

#====================================================================================== -

#> pls ---- 
indx <- createMultiFolds(traindat$trt, k = 10, times = 10)
ctrl <- caret::trainControl(method = "repeatedcv", 
                            index = indx,
                            savePredictions = TRUE,
                            verboseIter = TRUE,
                            selectionFunction = "oneSE")
plsr <- caret::train(trt ~ .,
                     data = traindat,
                     preProcess = c("center", "scale"),
                     method = "pls",
                     tuneLength = 15, 
                     trControl = ctrl,
                     importance = TRUE)
saveRDS(plsr, paste0(sinkdir, "pls.rds"))

#> rf-ranger ----
mtry <- c(1, 2, 5, 9, 14, 20, 30, 45, 70, 100, 200)
min.node.size = c(1, 2, 5)
tune_grid <- expand.grid(mtry = mtry,
                         splitrule = "gini", #default
                         min.node.size = min.node.size)
ctrl <- caret::trainControl(method = "repeatedcv",
                            number = 10,
                            rep = 10,
                            classProbs = TRUE,
                            verbose = TRUE)
rf_ranger <- caret::train(trt ~ .,
                          data = traindat,
                          preProc = c("center", "scale"),
                          method = "ranger",
                          tuneGrid = tune_grid,
                          importance = "permutation",
                          num.trees = 2000,
                          trControl = ctrl)
saveRDS(rf_ranger, paste0(sinkdir, "rf.rds"))

#====================================================================================== -

#Model generalization ----

##Model performance is evaluated on plots of the FPWW022 experiments
##not used in model training

#load models 
mod_names <- as.list(list.files(sinkdir, pattern = ".rds"))
mods <- lapply(paste0(sinkdir, mod_names), readRDS)
names(mods) <- gsub(".rds", "", mod_names)

#create predictions and class probabilities
preds_prob <- lapply(mods, predict.train, newdata = testdat[ , names(testdat) != "trt"], type = "prob", na.action = na.pass)
preds_class <- lapply(mods, predict.train, newdata = testdat[ , names(testdat) != "trt"], type = "raw")

#plots included in the training set
train_plots <- data_FPWW023 %>% dplyr::filter(grepl("FPWW022", Plot_ID)) %>% pull(Plot_ID) %>% unique()
#extract plots of the main experiment not used for training
plot_seq <- data_FPWW022 %>% 
  dplyr::filter(!Plot_ID %in% train_plots) %>% 
  mutate(Gen_Name = ifelse(is.na(Gen_Name), "not known", Gen_Name)) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  pull(Plot_ID)
#get design
plot_inf <- data_FPWW022 %>% 
  mutate(trt = as.factor(trt)) %>% 
  dplyr::select(-contains("SI_")) %>%
  dplyr::filter(Plot_ID %in% plot_seq)

#Create full dataframe for both classification algorithms
preds <- list()
for(i in names(mods)){
  preds[[i]] <- cbind(plot_inf, preds_prob[i], preds_class[i]) %>% 
    tibble::as_tibble() %>% 
    rename(prob_ctrl = 11, prob_dis = 12, pred = 13) %>% 
    mutate(acc = ifelse(pred == "ctrl", "TRUE", "FALSE") %>% as.factor(),
           #add count if wrongly predicted as "diseased"
           count = ifelse(pred == "dis", 1, 0))
}
preds <- bind_rows(preds, .id = "model")

#Distribution of class probabilities
distr <- ggplot(preds) +
  geom_histogram(aes(prob_ctrl), bins = 50) +
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_continuous(limits = c(0, 1)) +
  geom_vline(xintercept = 0.5, col = "red") +
  facet_wrap(~model) +
  theme_bw()

#Spatial distribution of class probablities
preds <- split(preds, preds$model)
desplot <- list()
for(i in names(preds)){
  desplot[[i]] <- desplot(prob_ctrl ~ RangeLot + RowLot, 
        data = preds[[i]], cex = 0.6, ticks = TRUE, show.key = TRUE,
        midpoint = 0.5,
        col.regions = RedGrayBlue,
        main = FALSE)
}

pdf(paste0(sinkdir, "class_dyn_distr_spat.pdf"), 8.5, 8.5)
gridExtra::grid.arrange(distr, desplot[[1]], desplot[[2]], ncol = 2)
dev.off()

#Accuarcy of models
tpr <- list()
for (i in names(mods)){
  true <- length(preds_class[[i]][preds_class[[i]] == "ctrl"])
  false <- length(preds_class[[i]][!preds_class[[i]] == "ctrl"])
  tpr[[i]] <- true/(true+false)
}

#====================================================================================== -
#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# FIT (S)PLS-DA MODELS FOR CLASSIFICATION INTO DISEASED / DISEASE-FREE PLOTS
# EXTRACT AND PLOT VARIABLE IMPORTANCE
# ASSESS "INTERNAL" MODEL PERFORMANCE VIA CV
# ASSESS "GENERAL" MODEL PERFORMANCE ON FPWW022 DATA
# ASSESS "DYNAMIC" MODEL PERFORMANCE ON FPWW022 DATA OF OTHER MEASUREMENT DATES

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
library(prospectr)
library(stringr)
#mixOmics has been removed from CRAN
#This script can be run exclusively with version 6.3.2. of mixOmics;
#(functions 'predict' and 'vip' will not work otherwise)
#This version of the package needs to be installed from the CRAN archive
#If these two function cause problems, 
#uninstalling and re-installing the package normally solves them
require(devtools)
install_version("mixOmics", version = "6.3.2", repos = "http://cran.us.r-project.org", 
                upgrade = "never", INSTALL_opts = c("--no-lock"))
library(mixOmics)

path_to_utils <- "" #Specify path to utils folder
path_to_data <- "" #Specify path to data folder

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT/plsda"))){
  dir.create(paste0(path_to_data, "OUT/plsda"), recursive = TRUE)
}
sinkdir <- paste0(path_to_data, "OUT/plsda/") #Specify output directory

source(paste0(path_to_utils, "002_pls_utils.R"))

#Load data
spc <- readRDS(paste0(path_to_data, "prep_data/spectra_plsda_bin6.rds")) %>% split(., .$meas_date)

#====================================================================================== -

#Fit models ----

if(!dir.exists(paste0(path_to_data, "OUT/plsda/Function_output"))){
  dir.create(paste0(path_to_data, "OUT/plsda/Function_output"))
}

result <- lapply(spc, f_splsda, sinkdir = sinkdir)
saveRDS(result, paste0(sinkdir, "Function_output/err_impvars.rds"))

#====================================================================================== -

#Load models ----

#get model names
mod_id <- as.list(list.files(paste0(sinkdir, "Models"), pattern = "plsda", recursive = TRUE)) %>% 
  lapply(., function(x) gsub(".rds", "", x))

#load models
mod <- as.list(list.files(paste0(sinkdir, "Models"), pattern = "plsda", recursive = TRUE, full.names = TRUE)) %>%
  lapply(readRDS)

#name models
names(mod) <- mod_id

#get dates
meas_dates <- sapply(mod_id, function(x) unlist(stringr::str_split(x, "_"))[2])

#get model info
mod_inf <- sapply(mod_id, function(x) unlist(stringr::str_split(x, "_"))[1])

#reorder list according to the dates
dates <- list()
for(i in meas_dates){
  dates[[i]] <- mod[grep(i, names(mod))]
}

#====================================================================================== -

#Create vip plots ----
#calculate vip
vip <- lapply(mod, mixOmics::vip)
#transform to long data frame for plotting
d <-lapply(vip, data.frame) %>%
  lapply(rownames_to_column) %>% 
  lapply(as.tibble) %>% 
  lapply("[", 1:4) %>% 
  purrr::map2(names(.), ~ add_column(.x, mod.id = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>% 
  rowwise() %>% 
  mutate(mod = grep("pls", unlist(str_split(mod.id, "_")), value = TRUE),
         meas_date = grep("2018", unlist(str_split(mod.id, "_")), value = TRUE),
         wvlt = gsub("rflt_", "", rowname) %>% as.numeric()) %>% 
  ungroup() %>%
  dplyr::select(wvlt, contains("comp"), mod, meas_date) %>% 
  tidyr::gather(comp, vip, contains("comp")) %>% 
  add_spc_range()
#extract one type of model
d <- filter(d, mod == "plsda" & !meas_date %in% c("2018-05-30", "2018-06-19", "2018-07-07", "2018-07-09", "2018-07-12"))

#create Figure
if(!dir.exists(paste0(sinkdir, "Figures"))){
  dir.create(paste0(sinkdir, "Figures"))
}

pdf(file = paste0(sinkdir, "Figures/vip_plsda_allcomp.pdf"))
ggplot(d) +
  geom_line(aes(x = wvlt, y = vip, group = interaction(comp, spc_range) , color = comp), size = 0.5) + 
  geom_hline(yintercept = 0.8) +
  xlab("wavelength") + ylab("VIP") +
  facet_wrap(~meas_date) +
  theme_bw() +
  theme(legend.title = element_blank())
dev.off()

#====================================================================================== -

#Exctract ncomp ----

ncomp <- lapply(mod, "[[", "ncomp") %>% unlist(.) %>% t(.) %>% t(.) %>% 
  data.frame() %>% rownames_to_column() %>% 
  tibble::as_tibble() %>% 
  rename(ncomp = 2, mod = 1) %>% 
  mutate(modn = lapply(strsplit(mod, "_"), "[[", 1) %>% unlist(),
         train = lapply(strsplit(mod, "_"), "[[", 2) %>% unlist()) %>% 
  dplyr::select(modn, train, ncomp) %>% 
  rename(mod = 1) %>% 
  filter(mod == "plsda" & train %in% unique(d$meas_date))

#====================================================================================== -

#"Internal" Model performance ----

data <- readRDS(paste0(sinkdir, "/Function_output/err_impvars.rds"))

#extract performance data
perf_data <- lapply(data, "[[", 1)

get_err_rate <- function(d) {
  err <- d$overall[nrow(d$overall)]
}

preds <- lapply(perf_data, function(y) lapply(y, get_err_rate))
dates0 <- names(preds)

out <- list()
for(i in 1:length(preds)) {
  out[[i]] <- preds[[i]] %>% unlist(.) %>% t(.) %>% t(.) %>% 
    data.frame() %>% rownames_to_column() %>% 
    tibble::as_tibble() %>% 
    rename(acc_int = 2, mod = rowname) %>% 
    mutate(train = dates0[i],
           test = dates0[i],
           acc_int = 1-acc_int) %>% 
    dplyr::select(mod, train, test, acc_int)
}
out <- do.call("rbind", out)

acc_int <- out %>% filter(mod == "plsda" & train %in% unique(d$meas_date))

#====================================================================================== -

#"General" Model performance ----

#Load train and test data
dat_train <- readRDS(paste0(path_to_data, "prep_data/spectra_plsda_bin6.rds"))
dat_test <- readRDS(paste0(path_to_data, "prep_data/spectra_plsda_FPWW022_bin6.rds"))

dat <- split(dat_train, dat_train$meas_date)

#for each date: get model, 
#get data for the most close-by measurement of the main experiment,
#use model to predict main experiment class
#compute accuracy
acc <- list()
for(i in names(dates)) {
  #get closest date for the FPWW022 measurements
  closest_date <- dat_test$meas_date[which.min(abs(difftime(i, unique(dat_test$meas_date))))]
  #get plots used for training 
  data <- dat[[i]] 
  trainplots <- unique(data$Plot_ID) #get plots used for training
  #predictors used in training
  preds <- names(data)[grepl("rflt_", names(data))]
  #extract relevant data, exclude data used for training
  test_data <- dat_test %>% 
    #exclude plots used for training; extract for relevant date
    filter(!Plot_ID %in% trainplots & meas_date == closest_date) %>%
    #select only preds used in training
    dplyr::select(one_of(preds)) %>% as.matrix()
  
  #get corresponding models
  modx <- dates[[i]]
  #extract model performance
  out <- sapply(modx, assess_mod_perf, test_data) %>% 
    #tidy up
    data.frame() %>% 
    rownames_to_column() %>% rowwise() %>% tibble::as.tibble() %>% 
    mutate(mod = strsplit(rowname, "_") %>% lapply("[[", 1) %>% unlist(),
           train = strsplit(rowname, "_") %>% lapply("[[", 2) %>% unlist(),
           test = as.character(closest_date)) %>% 
    rename(acc_gen = 2) %>% 
    dplyr::select(mod, train, test, acc_gen)
  acc[[i]] <- out
}

perf <- do.call("rbind", acc)

acc_gen <- perf %>% filter(mod == "plsda" & train %in% unique(d$meas_date))

#to create unique factor level for plotting (see below)
acc_gen[9, 4] <- 0.7972223

#====================================================================================== -

#Create full VIP plots ----

#combine all model characteristics
all_data <- left_join(d, acc_int, by = c("mod" ,"meas_date" = "train")) %>% 
  full_join(., acc_gen, by = c("mod", "meas_date" = "train")) %>% 
  full_join(., ncomp, by = c("mod", "meas_date" = "train"))

#Create VIP plot including model characteristics and performance metric

pdf(paste0(sinkdir, "Figures/full_vip_plsda.pdf"), width = 8.5, height = 7)
ggplot(all_data) +
  geom_rect(aes(xmin=680, xmax=750, ymin=-Inf, ymax=Inf), fill = "gray85") +
  geom_line(aes(x = wvlt, y = vip, group = interaction(comp, spc_range) , color = comp), size = 0.5) + 
  geom_hline(yintercept = 0.8) +
  xlab("wavelength") + ylab("VIP") +
  annotate("text", x = 2000, y = 1.9, label = paste0("acc[int]", "==", round(unique(all_data$acc_int), 2)),
           size = 3, parse = TRUE) +
  annotate("text", x = 2000, y = 1.65, label = paste0("acc[gen]", "==", round(unique(all_data$acc_gen), 2)),
           size = 3, parse = TRUE) +
  annotate("text", x = 1250, y = 1.9, label = paste0("ncomp", "==", all_data$ncomp[seq(1, 10260/4, by = 285)]),
           size = 3, parse = TRUE) +
  theme_bw() +
  facet_wrap(~meas_date) +
  theme(legend.title = element_blank(),
        strip.text = element_text(face = "bold"))
dev.off()

#====================================================================================== -

#"Dynamic" Model performance: external ----
# Only specificity can be assessed, but not sensitivity!

train_dates <- names(dates)

acc_dyn_ext <- list()
#iterate over dates;
#use each date as training,
#check performance on all other dates
for (i in train_dates){
  print(i)
  train_date <- i
  #list all dates except training date
  testdates <- unique(dat_test$meas_date)
  #get corresponding models
  modx <- dates[[train_date]]
  #for each date: get model, 
  #get data for the most close-by measurement of the main experiment,
  #use model to predict main experiment class
  #compute accuracy
  acc <- list()
  for(j in 1:length(testdates)){
    print(j)
    test_date <- testdates[j]
    #extract newdata
    test_data <- dat_test %>% 
      #extract date and drop incomplete cases
      dplyr::filter(meas_date == test_date) %>%
      #get predictors
      dplyr::select(contains("rflt_")) %>% 
      as.matrix()
    #extract performance
    out <- sapply(modx, assess_mod_perf, test_data) %>% 
      #tidy up
      data.frame() %>% 
      rownames_to_column() %>% rowwise() %>% tibble::as.tibble() %>% 
      mutate(model = strsplit(rowname, "_") %>% lapply("[[", 1) %>% unlist(),
             train = strsplit(rowname, "_") %>% lapply("[[", 2) %>% unlist(),
             test = test_date) %>% 
      rename(perf = 2) %>% 
      dplyr::select(model, train, test, perf)
    acc[[j]] <- out
  }
  acc_dyn_ext[[i]] <- do.call("rbind", acc)
}

acc_dyn_ext <- do.call("rbind", acc_dyn_ext)

perf_dyn_ext <- acc_dyn_ext %>% 
  filter(model == "plsda") %>% 
  mutate(train = as.Date(train)) %>% 
  rowwise() %>% 
  filter(test >= train)

mdates <- unique(perf_dyn_ext$test)

ggplot(perf_dyn_ext) +
  geom_point(aes(x = test, y = perf, group = train, col = as.factor(train))) +
  geom_line(aes(x = test, y = perf, group = train, col = as.factor(train))) +
  geom_abline(intercept = 0.5, slope = 0, size = 1, linetype = "dashed") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_x_date(breaks = mdates,
               labels = mdates) +
  xlab("Test date") + ylab("Overall accuracy") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 7.5, angle = 45, vjust = 1, hjust = 1))

#====================================================================================== -

#"Dynamic" Model performance: internal ----

train_dates <- names(dates)

acc_dyn_int <- list()
#iterate over dates;
#use each date as training,
#check performance on all other dates
for (i in train_dates){
  print(i)
  train_date <- i
  #list all dates in the test set 
  #Note: here, the "training" set is used as a test set, but on other dates!
  testdates <- unique(dat_train$meas_date)
  #get corresponding models
  modx <- dates[[train_date]]
  #for each date: get model, 
  #get data for the most close-by measurement of the main experiment,
  #use model to predict main experiment class
  #compute accuracy
  acc <- list()
  for(j in 1:length(testdates)){
    test_date <- testdates[j]
    #extract newdata
    test_data <- dat_train %>% 
      #extract date and drop incomplete cases
      dplyr::filter(meas_date == test_date) %>%
      #get predictors
      dplyr::select(contains("rflt_")) %>% 
      as.matrix()
    #create vector of true class labels
    labels <- c(rep("ctrl", 36), rep("dis", 36))
    #one missing ctrl plot on 2018-06-20
    if (j == 7){
      labels <- c(rep("ctrl", 35), rep("dis", 36))
    }
    #ony one replicate measured on 2018-06-10
    if (j == 4){
      labels <- c(rep("ctrl", 18), rep("dis", 18))
    }
    #calculate model performance
    out <- sapply(modx, assess_mod_perf2, testdata = test_data, labels = labels) %>% 
      #tidy up
      data.frame() %>% 
      rownames_to_column() %>% rowwise() %>% tibble::as.tibble() %>% 
      mutate(model = strsplit(rowname, "_") %>% lapply("[[", 1) %>% unlist(),
             train = strsplit(rowname, "_") %>% lapply("[[", 2) %>% unlist(),
             test = test_date) %>% 
      rename(perf = 2) %>% 
      dplyr::select(model, train, test, perf)
    acc[[j]] <- out
  }
  acc_dyn_int[[i]] <- do.call("rbind", acc)
}

acc_dyn_int <- do.call("rbind", acc_dyn_int)

perf_dyn_int <- acc_dyn_int %>% 
  filter(model == "plsda") %>% 
  mutate(train = as.Date(train)) %>% 
  filter(!train %in% c(as.Date("2018-07-07"), as.Date("2018-07-09"), as.Date("2018-07-12"))) %>% 
  rowwise() %>% 
  filter(test >= train)

mdates <- unique(perf_dyn_int$test)

pdf(paste0(sinkdir, "Figures/dyn_perf.pdf"), width = 9, height = 4.5)
ggplot(perf_dyn_int) +
  geom_point(aes(x = test, y = perf, group = train, col = as.factor(train)), size = 2) +
  geom_line(aes(x = test, y = perf, group = train, col = as.factor(train)), size = 1, alpha = 0.5) +
  geom_abline(intercept = 0.5, slope = 0, size = 1, linetype = "dashed") +
  scale_color_hue(l=60, c=40) +
  scale_x_date(breaks = mdates,
               labels = mdates) +
  xlab("Test date") + ylab("Overall accuracy") +
  theme_bw(base_size = 15) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 7.5, angle = 45, vjust = 1, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
dev.off()

#====================================================================================== -

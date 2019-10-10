
#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Perform recursive feature elimination with rf and cubist as base-learners

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

list.of.packages <- c("doParallel", "Rmpi", "caret", "Cubist", "tidyverse", "ranger")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos='https://stat.ethz.ch/CRAN/')

.libPaths("")

library(tidyverse)
library(caret)
library(ranger)
library(Cubist)
library(doParallel)

path_to_data <- ""
path_to_utils <- ""

source(paste0(path_to_utils, "005_rfe_utils.R"))

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT/rfe/regr"))){
  dir.create(paste0(path_to_data, "OUT/rfe/regr"))
}
sinkdir <- paste0(path_to_data, "OUT/rfe/regr/") #Specify output directory

#====================================================================================== -

data <- readRDS(paste0(path_to_data, "OUT/SI/preds/predmat_regr_allctrl.rds"))

# The candidate set of the number of predictors to evaluate
subsets <- c(length(data), 200, 150, 120, 105, 90, 75, 
             60, 50, 40, 35, 30, 25, 20, 17, 14, 12:1)

# Set up a cluster
p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
print(paste0("Detected cores: ", p)) #check cluster set up
cluster <- makeCluster(p - 1, type = "MPI", outfile = "")
registerDoParallel(cluster)
clusterEvalQ(cluster, {
  library(MASS)
  library(ranger)
  library(Cubist)
  library(caret)
  library(tidyverse)
})

rfe <- perform_rfe(response = "sev", base_learner = "ranger", type = "regression",
                   p = 0.7, times = 3, groups = 9, 
                   subsets = subsets, data = data,
                   importance = "permutation",
                   num.trees = 1000)

stopCluster(cluster)
registerDoSEQ()

saveRDS(paste0(sinkdir, "rfe_regr_ranger.rds"))

#====================================================================================== -

out <- tidy_rfe_output(rfe, "ranger")

PROF <- plot_perf_profile(out)

#====================================================================================== -

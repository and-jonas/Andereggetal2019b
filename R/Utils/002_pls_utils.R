#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

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

#classification using (S)PLS-DA
f_splsda <- function(data, nrepeat = 50, sinkdir) {
  
  #iter
  print(paste("fitting ", unique(data$meas_date), sep = ""))
  
  #get data per measurement date
  sub <- data %>% 
    dplyr::select(trt, contains("rflt"))

  #predictor matrix
  X <- sub %>% 
    dplyr::select(starts_with("rflt_")) %>%
    as.matrix()
  
  #response vector
  Y <- sub %>%
    dplyr::select(trt) %>%
    as.matrix() %>%
    as.factor()
  
  ## 1) PLS-DA; determine ncomp ----
  
  #fit plsda model with a maximum of 20 components
  plsda <- mixOmics::plsda(X, Y, ncomp = 20)
  
  #evaluate the plsda model using 10 fold cross-validation, repeated 10 times
  #and identify the optimal number of components to choose
  perf.plsda <- mixOmics::perf(plsda, validation = "Mfold", folds = 10, 
                               progressBar = FALSE, auc = TRUE, nrepeat = 10)
  
  #create performance plot
  plot(perf.plsda)
  plot01 <- recordPlot()
  plot.new() ## clean up device
  
  #extract the absolut best tuning parameter setting
  min.err <- min(perf.plsda$error.rate$overall[,"max.dist"])
  bestTune <- as.numeric(which.min(perf.plsda$error.rate$overall[,"max.dist"]))
  
  #extract the corresponding sd and calculate threshold
  sd.min.err <- perf.plsda$error.rate.sd$overall[bestTune,"max.dist"]
  threshold <- min.err + sd.min.err
  
  #extract lowest tuning parameter included in the tolerance margin of 1 sd
  #alternatively, the $choice.ncomp could be used 
  #(optimal number based on t-tests that test for a significant difference in the mean error rate between components)
  optTune <- as.numeric(min(which(perf.plsda$error.rate$overall[,"max.dist"] <= threshold)))

  ## 2) Final PLS-DA model; assess performance ----
  
  plsda <- mixOmics::plsda(X, Y, ncomp = optTune)
  #save model
  if(!exists(paste0(sinkdir, "Models"))){
    dir.create(paste0(sinkdir, "Models"))
  }
  saveRDS(plsda, paste(sinkdir, "Models/plsda_", 
                       unique(data$meas_date), ".rds", sep = ""))
  
  perf <- perf(plsda, validation = "Mfold", folds = 5,
               dist = 'max.dist', nrepeat = nrepeat, 
               progressBar = FALSE) 
  
  #create performance plot
  plot(perf)
  plot02 <- recordPlot()
  plot.new() ## clean up device
  
  #create ROC curve
  tail(auroc(plsda, roc.comp = optTune), 1)
  plot03 <- recordPlot()
  plot.new() ## clean up device
  
  #get error rates
  err_rate_plsda <- perf$error.rate
  
  #=============================================================================================================== -
  
  ## 3) SPLSDA. Tune number of variables to use in each component (done sequentially) ----
  
  # grid of possible keepX values that will be tested for each comp 
  list.keepX <- c(1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 75)
  
  #based on the cross validation results above, we set ncomp = optTune
  tune.splsda <- tune.splsda(X, Y, ncomp = optTune, validation = 'Mfold', folds = 5, 
                             progressBar = FALSE, dist = 'max.dist', 
                             test.keepX = list.keepX, nrepeat = nrepeat) #nrepeat 50-100 for better estimate
  
  #create performance plot
  plot(tune.splsda)
  plot10 <- recordPlot()
  plot.new() ## clean up device
  
  error <- tune.splsda$error.rate  # error rate per component for the keepX grid
  ncomp <- tune.splsda$choice.ncomp$ncomp # optimal number of components based on t-tests
  
  select.keepX <- tune.splsda$choice.keepX[1:ncomp]  # optimal number of variables to select
  
  ## 4) Fit final SPLSDA model and assess performance ----
  
  #fit the final model using only 1 variable for each dimension
  splsda <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX) 
  #save model
  saveRDS(splsda, paste(sinkdir, "Models/splsda.full_", 
                        unique(data$meas_date), ".rds", sep = ""))
  
  #create ROC curve
  tail(auroc(splsda, roc.comp = ncomp), 1)
  plot11 <- recordPlot()
  plot.new() ## clean up device
  
  #the performance of the final splsda model is assessed with the perf function 
  #by specifying a prediction distance
  perf <- perf(splsda, validation = "Mfold", folds = 5,
               dist = 'max.dist', nrepeat = nrepeat,
               progressBar = FALSE) 
  
  #create performance plot
  plot(perf) 
  plot12 <- recordPlot()
  plot.new() ## clean up device
  
  #extract stable features
  imp_vars_all <- perf$features$stable
  
  #get error rates
  err_rate_splsda_full <- perf$error.rate
  
  #=============================================================================================================== -
  
  ## 5) Fit SPLSDA model, force number variables for each component to one and assess performance ----
  
  #force number of variables used per component to 1
  select.keepX[select.keepX < 200] <- 1
  
  #fit the final model using only 1 variable for each dimension
  splsda <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX) 
  #save model
  saveRDS(splsda, paste(sinkdir, "Models/splsda.red_", 
                        unique(data$meas_date), ".rds", sep = ""))
  
  #create ROC curve
  tail(auroc(splsda, roc.comp = ncomp), 1)
  plot20 <- recordPlot()
  plot.new() ## clean up device
  
  #Assess performance of the final splsda model
  perf <- perf(splsda, validation = "Mfold", folds = 5,
               dist = 'max.dist', nrepeat = nrepeat,
               progressBar = FALSE) 
  
  #create performance plot
  plot(perf) 
  plot21 <- recordPlot()
  plot.new() ## clean up device
  
  imp_vars_red <- perf$features$stable
  
  #get error rates
  err_rate_splsda_red <- perf$error.rate
  
  #=============================================================================================================== -  
  
  ## 6) CREATE FINAL FUNCTION OUTPUT ----
  
  #generate final graph
  if(!exists(paste0(sinkdir, "Plots"))){
    dir.create(paste0(sinkdir, "Plots"))
  }
  pdf(paste(sinkdir, "Plots/results_", unique(data$meas_date), ".pdf", sep = ""), onefile=TRUE)    
  replayPlot(plot01)
  replayPlot(plot02)
  replayPlot(plot03)
  replayPlot(plot10)
  replayPlot(plot11)
  replayPlot(plot12)
  replayPlot(plot20)
  replayPlot(plot21)
  
  dev.off()
  
  #generate list of outputs: 
  #error rates of the different models
  err_rates <- list(err_rate_plsda, err_rate_splsda_full, err_rate_splsda_red)
  names(err_rates) <- c("plsda", "splsda_full", "splsda_red")
  #important variables
  imp_vars <- list(imp_vars_all, imp_vars_red)
  names(imp_vars) <- c("full", "red")
  #final function output
  out <- list(err_rates, imp_vars)

  return(out)
  
}

#Assess model performance on test data;
#data of the main experiment, which are all "healthy"
assess_mod_perf <- function(model, testdata){
  #create predictions for test seet
  test.predict <- predict(model, testdata)
  Prediction <- test.predict$class$max.dist[, ncol(test.predict$class$max.dist)]
  #compute accuracy of predictions
  tn <- length(Prediction[Prediction == "ctrl"]) #true negatives
  fp <- length(Prediction[Prediction == "dis"]) #false positives
  acc <- tn/nrow(testdata)
  return(acc)
}

#Assess model performance on test data; 
#data from the FPWW023 experiment
assess_mod_perf2 <- function(model, testdata, labels){
  ##model: a (s)plsda object from mixOmics
  ##testdata: a matrix containing the spectral data of plots to be classified
  ##labels: a vector holding the true class labels in the same order as the test data
  #create predictions for test set
  test.predict <- predict(model, testdata)
  #Extract prediction from the model with the optimal number of components
  Prediction <- test.predict$class$max.dist[, ncol(test.predict$class$max.dist)]
  #compute accuracy of predictions
  predobs <- cbind(labels, Prediction) %>% as.data.frame(.)
  true <- nrow(predobs[as.character(predobs$labels) == as.character(predobs$Prediction),])
  acc <- true/nrow(testdata)
  return(acc)
}

#====================================================================================== -
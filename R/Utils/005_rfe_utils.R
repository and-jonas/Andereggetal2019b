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

#Perform recursive feature elimination
perform_rfe <- function(response, base_learner = "ranger", type = "regression",
                        p = 0.75, times = 30, groups = 9, 
                        subsets, data,
                        ...) {
  
  #create multifolds for repeated n-fold cross validation
  index <- caret::createDataPartition(pull(data[response]), p = p, times = times, groups = ifelse(is.numeric(groups), groups, 2))

  #outer resampling
  #CV of feature selection
  out <- list()
  for(i in 1:length(index)){
    
    #Verbose
    print(paste("resample ", i, "/", length(index), sep = ""))
    
    #use indices to create train and test data sets for the resample
    ind <- as.numeric(index[[i]])
    train <- data[ind,]
    test <- data[-ind, ]
    
    #for each subset of decreasing size
    #tune/train rf and select variables to retain
    keep_vars <- drop_vars <- test_perf <- train_perf <- npred <- NULL
    for(j in 1:length(subsets)){
      
      #define new training data
      #except for first iteration, where the full data set ist used
      if(exists("newtrain")) {train = newtrain}
      
      #Verbose iter
      print(paste("==> subset size = ", length(train)-1, sep = ""))
      
      #define tune grid
      if(base_learner == "ranger"){
        #adjust mtry parameter to decreasing predictor set
        #maximum mtry at 200
        mtry <- ceiling(seq(1, length(train[-1]), len = 7)) %>% unique()
        if(any(mtry > 250)){
          mtry <- mtry[-which(mtry >= 250)]
        }
        min.node.size <- c(5)
        tune_grid <- expand.grid(mtry = mtry,
                                 splitrule = ifelse(type == "regression", "variance", "gini"),
                                 min.node.size = ifelse(type == "regression", 5, 1)) 
      } else if(base_learner == "cubist"){
        tune_grid <- expand.grid(committees = c(1, 2, 5, 10),
                                 neighbors = c(0))
      }
      
      #define inner resampling procedure
      ctrl <- caret::trainControl(method = "repeatedcv",
                                  number = 10,
                                  rep = 1,
                                  verbose = FALSE,
                                  allowParallel = TRUE,
                                  savePredictions = TRUE,
                                  classProbs = ifelse(type == "classification", TRUE, FALSE))
      
      #define model to fit
      formula <- as.formula(paste(response, " ~ .", sep = ""))
      
      #tune/train random forest
      fit <- caret::train(formula,
                          data = train,
                          preProc = c("center", "scale"),
                          method = base_learner,
                          tuneGrid = tune_grid,
                          trControl = ctrl,
                          ...)
      
      if(type == "regression"){
        #extract predobs of each cv fold
        predobs_cv <- plyr::match_df(fit$pred, fit$bestTune, on = names(fit$bestTune))
        #Average predictions of the held out samples;
        predobs <- predobs_cv %>% 
          group_by(rowIndex) %>% 
          dplyr::summarize(obs = mean(obs),
                           mean_pred = mean(pred))
        #get train performance
        train_perf[j] <- caret::getTrainPerf(fit)$TrainRMSE
        #get test performance
        test_perf[j] <- rmse(test %>% pull(response), caret::predict.train(fit, test))
      } else if (type == "classification"){
        #get train accuracy
        train_perf[j] <- caret::getTrainPerf(fit)$TrainAccuracy
        #get test accuracy
        test_perf[j] <- get_acc(fit, test)
      }
      
      #number of preds used
      npred[[j]] <- length(train)-1
      
      #extract retained variables
      #assign ranks
      #define reduced training data set
      if(j < length(subsets)){
        #extract top variables to keep for next iteration
        keep_vars[[j]] <- varImp(fit)$importance %>% 
          tibble::rownames_to_column() %>% 
          as_tibble() %>% dplyr::rename(var = rowname) %>%
          arrange(desc(Overall)) %>% slice(1:subsets[j+1]) %>% pull(var)
        #extract variables dropped from dataset
        drop_vars[[j]] <- names(train)[!names(train) %in% c(keep_vars[[j]], response)] %>% 
          tibble::enframe() %>% mutate(rank = length(subsets)-j+1) %>% 
          dplyr::select(value, rank) %>% dplyr::rename(var = value)
        #define new training data
        newtrain <- dplyr::select(train, response, keep_vars[[j]])
        #last iteration
      } else {
        drop_vars[[j]] <- names(train)[names(train) != response] %>% 
          tibble::enframe() %>% mutate(rank = length(subsets)-j+1) %>% 
          dplyr::select(value, rank) %>% rename(var = value)
      }
    } #END OF FEATURE ELIMINATION ON RESAMPLE i
    #clean environment 
    rm("newtrain")
    #gather results for resample i
    ranks <- drop_vars %>% do.call("rbind", .)
    out[[i]] <- list(ranks, train_perf, test_perf, npred)
  } #END OF OUTER RESAMPLING
  return(out)
}

#Create a tidy output
tidy_rfe_output <- function(data, base_learner){
  #tidy up list output
  subsets <- data[[1]][[length(data[[1]])]]
  ranks <- lapply(data, "[[", 1) %>% 
    Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = "var"), .) %>% 
    purrr::set_names(., c("var", paste("Resample", 1:length(data), sep = "")))
  RMSEtrain <- lapply(data, "[[", 2) %>% lapply(., cbind, subsets) %>%
    lapply(., as_tibble) %>% Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = "subsets"), .) %>% 
    dplyr::select(subsets, everything()) %>% 
    purrr::set_names(c("subset_size", paste("Resample", 1:length(data), sep = "")))
  RMSEtest <- lapply(data, "[[", 3) %>% lapply(., cbind, subsets) %>%
    lapply(., as_tibble) %>% Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = "subsets"), .) %>% 
    dplyr::select(subsets, everything()) %>% 
    purrr::set_names(c("subset_size", paste("Resample", 1:length(data), sep = "")))
  #average across resamples, get sd and means
  Trainperf <- RMSEtrain %>%
    gather(resample, RMSE, contains("Resample")) %>%
    group_by(subset_size) %>%
    arrange(subset_size) %>%
    summarise_at(vars(RMSE), funs(mean, sd), na.rm = TRUE) %>%
    mutate(set = "Train")
  Testperf <- RMSEtest %>% 
    gather(resample, RMSE, contains("Resample")) %>%
    group_by(subset_size) %>%
    arrange(subset_size) %>%
    summarise_at(vars(RMSE), funs(mean, sd), na.rm = TRUE) %>%
    mutate(set = "Test")
  Perf <- bind_rows(Trainperf, Testperf) %>% mutate(algorithm = base_learner)
  #average ranks
  robranks <- ranks %>% 
    gather(resample, rank, contains("Resample")) %>%
    group_by(var) %>%
    summarise_at(vars(rank), funs(mean, sd), na.rm = TRUE) %>%
    arrange(mean)
  tidy_out <- list(Perf, robranks)
  return(tidy_out)
}

#Plot performance profiles
plot_perf_profile <- function(data){
  pd <- position_dodge(0.5) # move them .05 to the left and right
  #plot performance profiles
  ggplot(data, aes(x = subset_size, y = mean, group = set, colour = set)) +
    geom_point(position = pd) + geom_line() +
    geom_errorbar(position = pd, aes(ymin = mean - sd, ymax = mean + sd), width = 1, alpha = 0.5) + 
    xlab("#Features") + ylab("RMSE") +
    scale_x_continuous(limits = c(-2.5, 32.5)) +
    facet_wrap(~algorithm) +
    theme_bw() +
    theme(legend.title=element_blank(),
          plot.title = element_text(size = 15, face = "bold"),
          strip.text = element_text(face = "bold"))
}

#====================================================================================== -

#Helper function to calculate accuracy
get_acc <- function(model, testdata) {
  preds_class <- caret::predict.train(model, newdata = testdata[ , names(testdata) != "trt"])
  true_class <- testdata$trt
  res <- cbind(preds_class, true_class) %>% data.frame()
  match <- ifelse(res$preds_class == res$true_class, 1, 0) %>% sum()
  acc <- match/nrow(testdata)
}

#Helper function to get RMSE
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted) ^ 2))
}

#====================================================================================== -
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

assess_mod_perf <- function(obj, plot_seq_inf){
  #record number of traning observations
  n_train <- nrow(obj$trainingData)
  #extract the performance metrics from the optimal model
  perf <- obj$resample %>% 
    summarize(mean_RMSE = mean(RMSE),
              sd_RMSE = sd(RMSE),
              mean_Rsq = mean(Rsquared),
              sd_Rsq = sd(Rsquared))
  #get the tuning parameters
  tune_pars <- names(obj$bestTune)
  bestTune <- obj$bestTune
  #extract predictions and corresponding 'real' observations
  #from the held out samples of each resampling iteration
  #for the optimal model
  predobs_cv <- plyr::match_df(obj$pred, obj$bestTune, on = tune_pars)
  #adjust plotseqinf if necessary
  if(n_train < 70){
    plot_seq_inf <- plot_seq_inf %>% dplyr::filter(Exp == "FPWW023")
  }
  #Average predictions of the held out samples;
  #add plot information and method
  predobs <- predobs_cv %>% 
    group_by(rowIndex) %>% 
    summarize(obs = mean(obs), 
              mean_pred = mean(pred),
              sd_pred = sd(pred),
              sem = sem_ci(pred)) %>% 
    bind_cols(., plot_seq_inf) %>% 
    mutate(algorithm  = obj$method)
  #extract pred and obs for further evaluations
  pred <- predobs$mean_pred
  obs <- predobs$obs
  # Squared bias
  SB = (mean(obs - pred, na.rm = TRUE))^2
  # Non-unity slope
  NU = mean((pred - mean(pred))^2) * (1 - lm(obs ~ pred)$coefficients[2])^2
  # Lack of correlation
  LC = mean((obs - mean(obs))^2) * (1 - cor(obs, pred, use = "pairwise.complete.obs")^2)
  # Proportional contributions of SB, NU and LC to MSE in percent
  SB_prop = round((mean(obs - pred, na.rm = TRUE))^2
                  / mean((pred - obs)^2) * 100, 0)
  NU_prop = round(mean((pred - mean(pred))^2)
                  * (1 - lm(obs ~ pred)$coefficients[2])^2 / mean((pred - obs)^2) * 100, 0)
  LC_prop = round(mean((obs - mean(obs))^2)
                  * (1 - cor(obs, pred, use = "pairwise.complete.obs")^2)
                  / mean((pred - obs)^2) * 100, 0)
  
  #======================================================================================== -
  
  #ADJUSTED
  #exclude predobs from ctrls
  predobs_red <- predobs %>% filter(grepl("FPWW023", Plot_ID)) #!!!! modified from filter(Exp == "FPWW023) !!!!
  #overwrite pred and obs
  pred <- predobs_red$mean_pred
  obs <- predobs_red$obs
  #calculate adjusted performance metrics
  RMSE_adj <- sqrt(sum((pred-obs)^2)/nrow(predobs_red))
  sd_RMSE_adj <- perf[,"sd_RMSE"]
  Rsq_adj <- summary(lm(obs ~ pred))$r.squared
  sd_Rsq_adj <- perf[,"sd_Rsq"]
  # Squared bias
  SB = (mean(obs - pred, na.rm = TRUE))^2
  # Non-unity slope
  NU = mean((pred - mean(pred))^2) * (1 - lm(obs ~ pred)$coefficients[2])^2
  # Lack of correlation
  LC = mean((obs - mean(obs))^2) * (1 - cor(obs, pred, use = "pairwise.complete.obs")^2)
  # Proportional contributions of SB, NU and LC to MSE in percent
  SB_prop_adj = round((mean(obs - pred, na.rm = TRUE))^2
                  / mean((pred - obs)^2) * 100, 0)
  NU_prop_adj = round(mean((pred - mean(pred))^2)
                  * (1 - lm(obs ~ pred)$coefficients[2])^2 / mean((pred - obs)^2) * 100, 0)
  LC_prop_adj = round(mean((obs - mean(obs))^2)
                  * (1 - cor(obs, pred, use = "pairwise.complete.obs")^2)
                  / mean((pred - obs)^2) * 100, 0)
  
  #======================================================================================== -
  
  #gather results
  out <- tibble("n_train" = n_train,
                "tune_pars" = list(tune_pars),
                "bestTune" = list(bestTune),
                "predobs" = list(predobs),
                "mean_RMSE" = perf[,"mean_RMSE"], "sd_RMSE" = perf[,"sd_RMSE"],
                "mean_Rsq" = perf[,"mean_Rsq"],"sd_Rsq" = perf[,"sd_Rsq"],
                "SB" = SB_prop,
                "NU" = NU_prop,
                "LC" = LC_prop,
                "predobs_adj" = list(predobs_red),
                "RMSE_adj" = RMSE_adj, "sd_RMSE_adj" = sd_RMSE_adj,
                "Rsq_adj"= Rsq_adj, "sd_Rsq_adj"= sd_Rsq_adj,
                "SB_adj" = SB_prop_adj,
                "NU_adj" = NU_prop_adj,
                "LC_adj" = LC_prop_adj)
  return(out)
}

plot_predobs <- function(perf, adjust = TRUE){
  #create labels
  predobs <- perf$predobs[[1]]
  label_rmse <- bquote(RMSE ==.(format(perf$mean_RMSE, digits = 2))*''%+-%''*.(format(perf$sd_RMSE, digits = 2)))
  label_rsq <- bquote(italic(R)^2 ==.(format(perf$mean_Rsq, digits = 2))*''%+-%''*.(format(perf$sd_Rsq, digits = 2)))
  label_trsamp <- paste0("n[train] == ", nrow(perf[["predobs"]][[1]]))
  #Adjust performance metrics if required
  if(adjust){
    #exclude predobs from ctrls
    predobs <- perf[["predobs_adj"]][[1]]
    #calculate adjusted performance metrics
    RMSE_adj <- sqrt(sum((predobs$mean_pred-predobs$obs)^2)/nrow(predobs))
    sd_RMSE_adj <- perf["sd_RMSE"] %>% pull()
    Rsq_adj <- summary(lm(predobs$obs ~ predobs$mean_pred))$r.squared
    sd_Rsq_adj <- perf["sd_Rsq"] %>% pull()
    #overwrite labels with corrected performance metrics
    label_rmse <- bquote(RMSE[adj]==.(format(RMSE_adj, digits = 2))*''%+-%''*.(format(sd_RMSE_adj, digits = 2)))
    label_rsq <- bquote(italic(R)[adj]^2 ==.(format(Rsq_adj, digits = 2))*''%+-%''*.(format(sd_Rsq_adj, digits = 2)))
  }
  label_tune <- vector()
  tune_pars <- perf$tune_pars[[1]]
  for(i in 1:length(tune_pars)){
    label_tune[i] <- paste0(tune_pars[i], " = ", as.character(round(as.numeric(perf[["bestTune"]][[1]][,tune_pars[i]]), 2)))
  }
  label_tune <- paste0(label_tune, collapse = "\n")
  #create the predobs plot
  cols <- c("1" = "red", "3" = "blue")
  predobs3 <- ggplot(predobs) +
    scale_color_manual(values = c("1" = "darkgreen", "3" = "darkblue")) +
    geom_point(aes(x = obs, y = mean_pred),
               shape = 16, size = 2, alpha = 0.5) +
    geom_errorbar(aes(x = obs, ymin = mean_pred - sem, ymax = mean_pred + sem), width = 0.005, alpha = 0.4) +
    geom_abline(intercept = 0, color = "red", size = 1, linetype="dashed") +
    xlab("Observed") + ylab("Predicted (repeated CV)") +
    geom_smooth(aes(x = obs, y = mean_pred), method = "lm") +
    scale_y_continuous(limits = c(0, 0.4)) +
    scale_shape_manual(values=c(1, 2)) +
    annotate("text", x = 0.1, y = 0.4, label = label_rmse,
             size = 3.5) +
    annotate("text", x = 0.1, y = 0.37, label = label_rsq,
             size = 3.5) +
    annotate("text", x = 0.1, y = 0.3225, label = label_trsamp,
             size = 3.25, parse = TRUE) +
    annotate("text", x = 0.1, y = 0.31, label = label_tune,
             size = 3.25, vjust = 1) +
    facet_wrap(~algorithm) +
    theme_bw(base_size = 15) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"))
}

#======================================================================================== -

# Function to calculate standard error of the mean
sem_ci <- function(x) {
  qt(0.975, df = length(na.omit(x)) - 1) *
    sqrt(var(x, na.rm = TRUE) / length(na.omit(x)))
}

#======================================================================================== -
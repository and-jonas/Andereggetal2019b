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

# wrapper function to extract model and dynamics parameters from parametric model fits
get_dyn_pars <- function(data, create_plot){
  dat <- data
  SI <- unique(dat$SI)
  #add model fits
  data_fits <- dat %>%
    as.data.frame() %>%
    tidyr::nest(-c(Plot_ID, SI)) %>%
    group_by(Plot_ID, SI) %>%
    mutate(fit_cgom = furrr::future_map(data, 
                                        ~ nls_multstart(value ~ Gompertz_constrained(b, M, tempsum = meas_GDDAH),
                                                        data = .x,
                                                        iter = 750,
                                                        start_lower = c(b = -0.2, M = 550),
                                                        start_upper = c(b = 0.1, M = 750),
                                                        convergence_count = 150, 
                                                        supp_errors = 'Y'))) %>% 
    tidyr::gather(met, fit, fit_cgom:fit_cgom)
  # new data frame of predictions
  new_preds <- dat %>%
    do(., data.frame(meas_GDDAH = seq(min(.$meas_GDDAH), max(.$meas_GDDAH), length.out = 100), stringsAsFactors = FALSE))
  # max and min for each curve
  max_min <- group_by(dat, Plot_ID) %>%
    summarise(., min_gGDDAH = min(meas_GDDAH), max_gGDDAH = max(meas_GDDAH)) %>%
    ungroup()
  # create new predictions
  preds2 <- data_fits %>%
    tidyr::unnest(fit %>% purrr::map(., broom::augment, newdata = new_preds)) %>%
    merge(., max_min, by = 'Plot_ID') %>%
    group_by(., Plot_ID, met) %>%
    filter(., meas_GDDAH > unique(min_gGDDAH) & meas_GDDAH < unique(max_gGDDAH)) %>%
    arrange(., Plot_ID, SI, meas_GDDAH) %>%
    rename(., value = .fitted) %>%
    ungroup()
  # check whether model converged
  convInfo <- data_fits %>%    
    #reconverto to wide, for now
    dplyr::filter(met == "fit_cgom") %>% 
    tidyr::spread(met, fit) %>%
    transmute(convInfo = purrr::map(purrr::map(fit_cgom, "convInfo"), "isConv"))
  #create plots
  if(create_plot == TRUE){
    #create sink directory
    if(!dir.exists(paste0(path_to_data, "OUT/SI/fits"))){
      dir.create(paste0(path_to_data, "OUT/SI/fits"))
    }
    ggplot() +
      geom_point(aes(meas_GDDAH, value), size = 0.3, shape = 1, dat) +
      geom_line(aes(meas_GDDAH, value, group = met, colour = met), alpha = 0.75, size = 0.3, preds2) +
      facet_wrap(~ Plot_ID, labeller = labeller(.multi_line = FALSE)) +
      scale_colour_manual(values = c("green4", "black", "red")) +
      scale_y_continuous(breaks = seq(0,10,2), limits = c(0, 10)) +
      geom_abline(slope = 0, intercept = 8.5) +
      geom_abline(slope = 0, intercept = 5.0) +
      geom_abline(slope = 0, intercept = 1.5) +
      theme_bw(base_size = 3.5, base_family = "Helvetica") +
      ylab("Index value") +
      xlab("°C days after heading") +
      ggtitle(unique(dat$SI))
    ggsave(paste0(sinkdir, "fits/", SI, ".pdf"),  width = 7, height = 5)
  }
  # get model parameters
  mod_pars <- data_fits %>%
    #reconverto to wide, for now
    dplyr::filter(met == "fit_cgom") %>% 
    tidyr::spread(met, fit) %>% 
    tidyr::unnest(fit_cgom %>% purrr::map(broom::tidy)) %>% 
    dplyr::select(1:4) %>% 
    tidyr::spread(term, estimate)
  # extract dynamics parameters from nls fits
  dynpars_gom <- preds2 %>% 
    #use only cgom for now
    dplyr::filter(met == "fit_cgom") %>% 
    group_by(Plot_ID, SI) %>%
    tidyr::nest() %>%
    mutate(pars_gom = purrr::map(data, extract_pars)) %>%
    dplyr::select(-data)
  # combine pars and convInfo
  mod_data <- list(dynpars_gom, convInfo, mod_pars) %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("Plot_ID", "SI")), .)
  # build output df containing all parameters and the design
  out_data <- mod_data %>% tidyr::unnest(pars_gom, convInfo) %>% 
    rowwise() %>%
    mutate(trt = ifelse(grepl("FPWW023", Plot_ID), "dis", "ctrl")) %>%
    ungroup() %>% 
    full_join(inf) %>%
    mutate(Lot = Lot %>% as.factor()) %>% 
    dplyr::select(Plot_ID, Exp, Lot, RangeLot, RowLot, Gen_Name, 
                  trt, convInfo, contains("heading"), everything()) %>% 
    tidyr::gather(par, val, b:dur2)
  return(out_data)
}

# Wrapper function around extract_pars to generate full parameter dataset
get_dyn_pars_lin <- function(data, create_plot = FALSE){
  dat <- data
  SI <- unique(dat$SI)
  #add model fits
  data_fits <- dat %>% as.data.frame() %>%
    tidyr::nest(-c(Plot_ID, SI)) %>%
    group_by(Plot_ID, SI) %>%
    mutate(fit_lin = purrr::map(data, lin_approx) %>% purrr::map(cbind.data.frame)) %>% 
    tidyr::gather(met, fit, fit_lin:fit_lin)
  preds2 <- data_fits %>%
    tidyr::unnest(fit) %>%
    group_by(., Plot_ID, met) %>%
    arrange(., Plot_ID, SI, meas_GDDAH) %>%
    rename(., value = .fitted) %>%
    ungroup()
  # extract dynamics parameters from nls fits
  dynpars_lin <- preds2 %>% 
    #use only cgom for now
    dplyr::filter(met == "fit_lin") %>% 
    group_by(Plot_ID, SI) %>%
    tidyr::nest() %>%
    mutate(pars_lin = purrr::map(data, extract_pars)) %>%
    dplyr::select(-data)
  # build output df containing all parameters and the design
  out_data <- dynpars_lin %>% tidyr::unnest(pars_lin) %>% 
    ungroup() %>% 
    full_join(inf, .) %>%
    mutate(Lot = Lot %>% as.factor()) %>% 
    tidyr::gather(par, val, t85:dur2) %>% 
    arrange(SI, Plot_ID)
  return(out_data)
}

# Extraction of parameters from fits
extract_pars <- function(data){
  t85 <- data[which(data[3] < 8.5)[1], "meas_GDDAH"]
  t50 <- data[which(data[3] < 5)[1], "meas_GDDAH"]
  t15 <- data[which(data[3] < 1.5)[1], "meas_GDDAH"]
  dur1 <- t50 - t85
  dur2 <- t15 - t85
  pars <- cbind(t85, t50, t15, dur1, dur2)
  names(pars) <- c("t85", "t50", "t15", "dur1", "dur2")
  return(pars)
}

# select and transform parameter data
select_pars <- function(data, which, pars){
  filter(data, par %in% pars) %>%
    rename(!!paste("SI", which, sep = "_") := SI) %>%
    mutate(par = paste(par, paste0(which, "ind"),  sep = "_")) %>%
    tidyr::spread(par, val) %>%
    dplyr::select(-one_of("convInfo"))
}

#====================================================================================== -

# Contrained Gompertz equation
Gompertz_constrained <- function(b, M, tempsum) {
  grenness_decay <- 10*exp(-exp(-b*(tempsum-M)))
  return(grenness_decay)
}

# Logistic equation
logistic <- function(A, C, b, M, tempsum) {
  grenness_decay <- A + C/(1+exp(-b*(tempsum-M)))
  return(grenness_decay)
}

# Flexible Gompertz equation
Gompertz_flex <- function(A, C, b, M, tempsum) {
  grenness_decay <- A + C*exp(-exp(-b*(tempsum-M)))
  return(grenness_decay)
}

# Linear interpolation
lin_approx <- function(data){
  out <- approx(data[,"meas_GDDAH"], data[,"value"], xout = seq(round(min(data[,"meas_GDDAH"], na.rm = TRUE), 0),
                                                                round(max(data[,"meas_GDDAH"], na.rm = TRUE), 0), 1))
  names(out) <- c("meas_GDDAH", ".fitted")
  out <- list(out)
  return(out)
}

#====================================================================================== -
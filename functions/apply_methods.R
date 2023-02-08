# functions to apply missing data methods

# set seed for both R and Python environments
# set.seed(19)
# reticulate::py_set_seed(19)

# complete case analysis
apply_CCA <- function(amp) {
  # list-wise deletion
  est <- na.omit(amp$amp) %>% 
    # fit regression 
    lm(Y ~ X1 + X2 + X3 + X4, .) %>% 
    # clean results
    broom::tidy(conf.int = TRUE) %>% 
    # choose estimates
    select(term, estimate, conf.low, conf.high) %>% 
    # add method name and missingness
    cbind(method = "CCA", mech = amp$mech, prop = amp$prop, .) 
  # output
  return(est)
}

# MICE imputation
apply_MICE <- function(amp) {
  # imputation with MICE
  imp <- mice::mice(amp$amp, method = "norm", maxit = 5, m = 5, printFlag = TRUE) #five iterations and five datasets
  # fit regression on each imputation
  est <- with(imp, lm(Y ~ X1 + X2 + X3 + X4)) %>% 
    # pool results
    mice::pool() %>% 
    # clean results
    broom::tidy(conf.int = TRUE) %>% 
    # select estimates
    select(term, estimate, conf.low, conf.high) %>% 
    # add method name
    cbind(method = "MICE", mech = amp$mech, prop = amp$prop, .)
  # output
  return(est)
}

# miceRanger imputation
# executes exactly the same functions as miceforest in Python
# miceRanger is essentially miceforest translated to R
apply_MICER <- function(amp) {
  # imputation with MICERanger
  mimp <- miceRanger::miceRanger(amp$amp, maxiter = 5, m = 5, verbose = TRUE) #five iterations and five datasets
  # generate m complete datasets
  comp <- miceRanger::completeData(mimp)
  # unlist lists from the complete() function add iteration and add incomplete data
  unlisted <- do.call(rbind, comp) %>% cbind(data.frame(.imp = rep(1:5 ,each = 200))) %>%
    rbind(cbind(.imp = 0, amp$amp), .)
  # convert to mids object, to apply Rubin's pooling rules
  imp <- as.mids(unlisted)
  # fit regression on each imputation
  est <- with(imp, lm(Y ~ X1 + X2 + X3 + X4)) %>% 
    # pool results
    mice::pool() %>% 
    # clean results
    broom::tidy(conf.int = TRUE) %>% 
    # select estimates
    select(term, estimate, conf.low, conf.high) %>% 
    # add method name
    cbind(method = "MICER", mech = amp$mech, prop = amp$prop, .)
  # output
  return(est)
}

# Python imputation

np <- import("numpy") #import numpy
pd <- import("pandas") #import pandas
sklearn_impute <- import("sklearn.impute") #import sklearn.impute
knn_imputer <- sklearn_impute$KNNImputer() #import KNNImputer from sklearn.impute as knn_imputer
it_imputer <- sklearn_impute$IterativeImputer() #import IterativeImputer from sklearn.impute as it_imputer

# KNNImputer imputataion
apply_KNN_imputer <- function(amp) {
  # KNNImputer
  knn_imputer <- sklearn_impute$KNNImputer() #default is 5 neighbours
  # coerce NaN's
  amp$amp[is.na(amp$amp)]<-NaN
  # impute data
  imp <- knn_imputer$fit_transform(amp$amp)
  dat <- as.data.frame(imp)
  # get estimates and rename columns because of Python syntax
  est <- rename(dat, Y = V1,  X1 = V2, X2 = V3, X3 = V4, X4 = V5) %>%
    lm(Y ~ X1 + X2 + X3 + X4, .) %>%
    # extract confidence intervals and select relevant estimates
    broom::tidy(conf.int = TRUE) %>%
    select(term, estimate, conf.low, conf.high) %>%
    cbind(method = "KNN", mech = amp$mech, prop = amp$prop, .)
  return(est)
}

# IterativeImputer imputataion
apply_IterativeImputer <- function(amp) {
  # IterativeImputer
  it_imputer <- sklearn_impute$IterativeImputer(max_iter = 5L, sample_posterior = TRUE,
                                              verbose = 2L) #default max iterations is 10
  # coerce NaN's                                            # using five iterations and five datasets
  amp_nan <- amp$amp
  amp_nan[is.na(amp_nan)]<-NaN
  # repeat m times
  repeat_times <- function(inc, m) {
    # impute data with IterativeImputer
    imp <- it_imputer$fit_transform(inc)
    dat <- as.data.frame(imp)
    # rename columns because of Python syntax
    dat <- rename(dat, Y = V1,  X1 = V2, X2 = V3, X3 = V4, X4 = V5)
    dat <- cbind(.imp = m, dat)
    return(dat)  
  }
  stacked_dfs <- purrr::map_dfr(1:5, ~repeat_times(amp_nan, m = .x)) %>%
    # add incomplete data with NA
    rbind(cbind(.imp = 0, amp$amp), .)
  # convert to a mids object
  imp <- as.mids(stacked_dfs)
  # get estimates and apply Rubin's pooling rules
  est <- with(imp, lm(Y ~ X1 + X2 + X3 + X4)) %>% 
    # pool results
    mice::pool() %>% 
    # clean results
    broom::tidy(conf.int = TRUE) %>% 
    # select estimates
    select(term, estimate, conf.low, conf.high) %>% 
    # add method name
    cbind(method = "ITIMP", mech = amp$mech, prop = amp$prop, .)
  # output
  return(est)
}

# MIDAS imputation function
apply_RMIDAS <- function(amp) {
  # imputation with MIDAS
  # preprocess data for MIDAS
  # amp$amp[is.na(amp$amp)]<-NaN
  conv <- rMIDAS::convert(amp$amp)
  # choose hyperparameters and impute with MIDAS algorithm
  # using the default tuning as described in the rMIDAS package/paper
  imp <- rMIDAS::train(conv, training_epochs = 10L, 
                       layer_structure = c(256, 256, 256),
                       input_drop = 0.8)
  # return m complete datasets
  comp <- rMIDAS::complete(imp, m = 5)
  # estimate regression models on m datasets with Rubin's rules
  r_est <- rMIDAS::combine("Y ~ X1 + X2 + X3 + X4", comp, family = stats::gaussian())
  # calculate confidence intervals
  conf.low <- r_est$estimate - qnorm(0.975) * r_est$std.error
  conf.high <- r_est$estimate + qnorm(0.975) * r_est$std.error
  # return estimates as a data frame
  est <- data.frame(method = "MIDAS", mech = amp$mech, prop = amp$prop, term = r_est$term, estimate = r_est$estimate, 
                    conf.low, conf.high)
  # output
  return(est)
}

# combine into one function
apply_methods <- function(amps, betas) {
  # apply CCA to each incomplete dataset
  CCA <- purrr::map_dfr(amps, ~{apply_CCA(.)})
  # impute with MICE and estimate effects
  MICE <-  purrr::map_dfr(amps, ~{apply_MICE(.)})
  # impute with miceRanger and estimate effects
  MICER <-  purrr::map_dfr(amps, ~{apply_MICER(.)})
  # impute with Python and estimate effects
  # KNNImputer
  KNN <-  purrr::map_dfr(amps, ~{apply_KNN_imputer(.)})
  # IterativeImputer
  ITIMP <-  purrr::map_dfr(amps, ~{apply_IterativeImputer(.)})
  # MIDAS
  MIDAS <-  purrr::map_dfr(amps, ~{apply_RMIDAS(.)})
  # combine estimates 
  ests <- rbind(CCA, MICE, MICER, KNN, ITIMP, MIDAS) %>%
    cbind(truth = c(0, betas))
  # output
  return(ests)
}

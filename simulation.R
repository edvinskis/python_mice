# Simulation script for ADS thesis
# Original script by Hanne Oberman

########################
### SETUP SIMULATION ###
########################

# R packages
library(dplyr)
library(mvtnorm)
library(mice)
library(miceadds)
# reticulate package for handling Python code in R
library(reticulate)
reticulate::use_condaenv("r-reticulate", required = TRUE) #virtual environment is user specific
# Python imputation (MIDAS imputation in R)
# for rMIDAS to work properly please use a Python environment that contains at least:
# python >3.7
# numpy > 1.21.6
# pandas >1.3.5
# scikit-learn >1.1.1
# tensorflow >2.7
# tensorflow-addons >0.17.0
# matplotlib > 3.5.2
library(rMIDAS)
# packages for plots
library(esquisse)
library(ggplot2)
library(ggh4x)
library(ggvis)

# functions
miceadds::source.all("./functions")

# randomness
# set seed for both R and Python environments
set.seed(19)
reticulate::py_set_seed(19)

# parameters
n_sim <- 1000
n_obs <- 200
betas <- c(-0.5, -0.1, 0.1, 0.5)
mis_mech = c("MCAR", "MAR")
mis_prop = c(0.1, 0.25, 0.5)

#################################
### TEST LOWER LEVEL FUCTIONS ###
#################################

# generate data
dat <- generate_complete(n_obs, betas)

# ampute data
amps <- induce_missingness(dat, mis_mech = "MCAR", mis_prop = 0.5)

# apply complete case analysis
CCA <- apply_CCA(amps[[1]])

# impute data with MICE
MICE <- apply_MICE(amps[[1]])

# impute data with miceRanger 
# executes exactly the same functions as miceforest in Python
# miceRanger is essentially miceforest translated to R
MICER <- apply_MICER(amps[[1]])

### impute data with Python ###
# impute data with KNNImputer
KNN <- apply_KNN_imputer(amps[[1]])

# impute data with IterativeImputer
ITIMP <- apply_IterativeImputer(amps[[1]])

# impute data with MIDAS
MIDAS <- apply_RMIDAS(amps[[1]])

##################################
### TEST HIGHER LEVEL FUCTIONS ###
##################################

amps <- create_data()
ests <- apply_methods(amps, betas)

################################
### COMBINE INTO ONE FUCTION ###
################################

simulate_once <- function(n_obs, betas, mis_mech, mis_prop) {
  # generate incomplete data
  amps <- create_data(
    sample_size = n_obs,
    effects = betas,
    mechanisms = mis_mech,
    proportions = mis_prop
  )
  # estimate regression coefficients
  ests <- apply_methods(amps, betas)
  # output
  return(ests)
}

################################
### TEST SIMULATION FUNCTION ###
################################

ests <- simulate_once(n_obs, betas, mis_mech, mis_prop)

######################
### RUN SIMULATION ###
######################

# repeat the simulation function n_sim times
results_raw <- replicate(
  n_sim, 
  simulate_once(n_obs, betas, mis_mech, mis_prop),
  simplify = FALSE
)
# save raw results
saveRDS(results_raw, "./results/raw_1000.rds")

########################
### EVALUATE RESULTS ###
########################
raw <- readRDS("./results/raw_1000.rds")

# calculate bias, coverage rate and CI width
performance <- evaluate_est(raw)

# simulation results across all conditions
perf_all_cond <- performance %>% 
  group_by(method) %>% 
  summarise(across(c(bias, cov, ciw), mean)) 

# simulation results split by condition
perf_by_cond <- performance %>% 
  group_by(method, mech, prop) %>% 
  summarise(across(c(bias, cov, ciw), mean))

# simulation results split by condition and regression coefficient
perf_by_cond_and_coeff <- performance %>% 
  group_by(method, mech, prop, term) %>% 
  summarise(across(c(bias, cov, ciw), mean)) 

# a rounded table for appendix
#table_appendix <- perf_by_cond_and_coeff %>% mutate(across(where(is.numeric), ~ round(., digits = 3)))
table_appendix <- perf_by_cond_and_coeff %>% mutate(across(where(is.numeric), ~ signif(., digits = 3)))
#readr::write_csv2(table_appendix,"performance_table.csv")

#############
### PLOTS ###
#############

##################################
### BIAS PLOT (X1, X2, X3, X4) ###
##################################

# 10% missingness
bias_plot_0.1 <- 
  performance %>%
  filter(!(term %in% "(Intercept)")) %>%
  filter(prop >= 0.1 & prop <= 0.1) %>%
  mutate(prop = recode(prop, "0.1" = "Missingness proportion = 0.10")) %>%
  mutate(mech = factor(mech, levels = c("MCAR", "MAR"))) %>%
  ggplot() +
  aes(x = bias, y = method, fill = mech, colour = term) +
  geom_vline(xintercept = 0, linetype="solid", color = "grey30") +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greys", direction = 1) +
  labs(
    # subtitle = "Bias of X1, X2, X3, X4 for different imputation methods 
    # under MCAR and MAR missingness mechanisms with a missingness proportion of 10%",
    x = "Bias",
    y = "Imputation method",
    fill = "Missingness mechanism",
    colour = "Term") +
  #coord_flip() +
  theme_bw() +
  theme(legend.position = "top") +
  facet_nested(vars(prop, term), scales = "free_x") + coord_cartesian(xlim = c(- 0.45, 0.5)) +
  theme(text=element_text(size = 11,  family = "mono"))

plot(bias_plot_0.1)

# 25% missingness
bias_plot_0.25 <- 
  performance %>%
  filter(!(term %in% "(Intercept)")) %>%
  filter(prop >= 0.25 & prop <= 0.25) %>%
  mutate(prop = recode(prop, "0.25" = "Missingness proportion = 0.25")) %>%
  mutate(mech = factor(mech, levels = c("MCAR", "MAR"))) %>%
  ggplot() +
  aes(x = bias, y = method, fill = mech, colour = term) +
  geom_vline(xintercept = 0, linetype="solid", color = "grey30") +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greys", direction = 1) +
  labs(
    # subtitle = "Bias of X1, X2, X3, X4 for different imputation methods
    # under MCAR and MAR missingness mechanisms with a missingness proportion of 25%",
    x = "Bias",
    y = "Imputation method",
    fill = "Missingness mechanism",
    colour = "Term") +
  #coord_flip() +
  theme_bw() +
  theme(legend.position = "top") +
  facet_nested(vars(prop, term), scales = "free_x") + coord_cartesian(xlim = c(- 0.45, 0.5)) +
  theme(text=element_text(size = 11,  family = "mono"))

plot(bias_plot_0.25)

# 50% missingness
bias_plot_0.5 <- 
  performance %>%
  filter(!(term %in% "(Intercept)")) %>%
  filter(prop >= 0.5 & prop <= 0.5) %>%
  mutate(prop = recode(prop, "0.5" = "Missingness proportion = 0.50")) %>%
  mutate(mech = factor(mech, levels = c("MCAR", "MAR"))) %>%
  ggplot() +
  aes(x = bias, y = method, fill = mech, colour = term) +
  geom_vline(xintercept = 0, linetype="solid", color = "grey30") +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greys", direction = 1) +
  labs(
    # subtitle = "Bias of X1, X2, X3, X4 for different imputation methods 
    # under MCAR and MAR missingness mechanisms with a missingness proportion of 50%",
    x = "Bias",
    y = "Imputation method",
    fill = "Missingness mechanism",
    colour = "Term") +
  #coord_flip() +
  theme_bw() +
  theme(legend.position = "top") +
  facet_nested(vars(prop, term), scales = "free_x") + coord_cartesian(xlim = c(- 0.45, 0.5)) +
  theme(text=element_text(size = 11,  family = "mono"))

plot(bias_plot_0.5)

######################################
### COVERAGE PLOT (X1, X2, X3, X4) ###
######################################

# Filtering for plots
cov_0.5 <- perf_by_cond_and_coeff %>%
  filter(!(term %in% "(Intercept)")) %>%
  filter(prop >= 0.5 & prop <= 0.5) %>%
  mutate(prop = recode(prop, "0.5" = "Missingness proportion = 0.50")) 
cov_0.25 <- perf_by_cond_and_coeff %>%
  filter(!(term %in% "(Intercept)")) %>%
  filter(prop >= 0.25 & prop <= 0.25) %>%
  mutate(prop = recode(prop, "0.25" = "Missingness proportion = 0.25")) 
cov_0.1 <- perf_by_cond_and_coeff %>%
  filter(!(term %in% "(Intercept)")) %>%
  filter(prop >= 0.1 & prop <= 0.1) %>%
  mutate(prop = recode(prop, "0.1" = "Missingness proportion = 0.10")) 

cov_for_plots <- rbind(cov_0.1,cov_0.25,cov_0.5)

cov_plots <- 
  cov_for_plots %>%
  mutate(mech = factor(mech, levels = c("MCAR", "MAR"))) %>%
  ggplot() +
  aes(x = method, y = cov, colour = mech) +
  geom_point(shape = 19) +
  scale_colour_manual(values = c("MCAR" = "gray74", "MAR" = "gray29")) +
  labs(
    # subtitle = 
    # "Coverage of X1, X2, X3, X4 for different imputation methods - 
    # under MCAR and MAR missingness mechanisms 
    # with three missingness proportions (10 %, 25% and 50%)",
    y = "Coverage",
    x = "Imputation method",
    colour = 
      "Missingness mechanism"
  ) +
  theme_bw() +
  theme(legend.position = "top") +
  facet_nested_wrap(vars(prop, term), scales = "fixed", ncol = 4) + 
  scale_y_continuous(breaks = seq(0.86, 0.97, by = 0.02)) +
  geom_hline(yintercept = 0.94, linetype="dashed", color = "darkgray") + 
  geom_hline(yintercept = 0.96, linetype="dashed", color = "darkgray") +
  geom_hline(yintercept = 0.95, linetype="dotted", color = "lightgray") +
  theme(text = element_text(size = 11,  family = "mono"))

plot(cov_plots)

#################################
### CIW PLOT (X1, X2, X3, X4) ###
#################################

ciw_0.5 <- performance %>%
  filter(!(term %in% "(Intercept)")) %>%
  filter(prop >= 0.5 & prop <= 0.5) %>%
  mutate(prop = recode(prop, "0.5" = "Missingness proportion = 0.50")) 
ciw_0.25 <- performance %>%
  filter(!(term %in% "(Intercept)")) %>%
  filter(prop >= 0.25 & prop <= 0.25) %>%
  mutate(prop = recode(prop, "0.25" = "Missingness proportion = 0.25")) 
ciw_0.1 <- performance %>%
  filter(!(term %in% "(Intercept)")) %>%
  filter(prop >= 0.1 & prop <= 0.1) %>%
  mutate(prop = recode(prop, "0.1" = "Missingness proportion = 0.10")) 

ciw_for_plots <- rbind(ciw_0.1, ciw_0.25, ciw_0.5)

ciw_plots <- 
  ciw_for_plots %>%
  mutate(mech = factor(mech, levels = c("MCAR", "MAR"))) %>%
  ggplot() +
  aes(x = ciw, y = method, fill = mech, colour = term) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greys", direction = 1) +
  labs(
    # subtitle = "Confidence interval width of X1, X2, X3, X4 for different imputation methods 
    # under MCAR and MAR missingness mechanisms with three missingness proportions",
    x = "Confidence interval width",
    y = "Imputation method",
    fill = "Missingness mechanism",
    colour = "Term") +
  theme_bw() +
  theme(legend.position = "top") +
  facet_nested_wrap(vars(prop, term), scales = "fixed", ncol = 4) + coord_cartesian(xlim = c(0.2, 1.06)) +
  theme(text = element_text(size = 11,  family = "mono"))

plot(ciw_plots)

###################
### PLOTS FOR X1###
###################

bias_plot_X1 <- 
  performance %>%
  filter(term %in% "X1") %>%
  #filter(mech %in% "MAR") %>% # if one would like to only look at one missingness mechanism
  filter(prop >= 0.5 & prop <= 0.5) %>%
  mutate(prop = recode(prop, "0.5" = "Missingness proportion = 0.5")) %>%
  mutate(mech = factor(mech, levels = c("MCAR", "MAR"))) %>%
  ggplot() +
  aes(y = bias, x = method, fill = mech, colour = term) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greys", direction = 1) +
  labs(
    # subtitle = "Bias of X1 for different imputation methods
    # under MCAR and MAR missingness mechanisms with a missingness proportion of 50%",
    y = "Bias",
    x = "Imputation method",
    fill = "Missingness mechanism",
    colour = "Term") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top") +
  facet_nested_wrap(vars(prop, term), scales = "fixed") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey30") +
  theme(text = element_text(size = 11,  family = "mono"))

plot(bias_plot_X1)

ciw_plot_X1 <- 
  performance %>%
  filter(term %in% "X1") %>%
  filter(prop >= 0.5 & prop <= 0.5) %>%
  mutate(prop = recode(prop, "0.5" = "Missingness proportion = 0.5")) %>%
  mutate(mech = factor(mech, levels = c("MCAR", "MAR"))) %>%
  ggplot() +
  coord_flip() +
  aes(x = ciw, y = method, fill = mech, colour = term) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Greys", direction = 1) +
  labs(
    # subtitle = "Confidence interval width of X1 for different imputation methods
    # under MCAR and MAR missingness mechanisms with a missingness proportion of 50%",
    x = "Confidence interval width",
    y = "Imputation method",
    fill = "Missingness mechanism",
    colour = "Term") +
  theme_bw() +
  theme(legend.position = "top") +
  facet_nested_wrap(vars(prop, term), scales = "fixed") +
  theme(text = element_text(size = 11,  family = "mono"))

plot(ciw_plot_X1)

cov_plot_X1 <- 
  perf_by_cond_and_coeff %>%
  filter(term %in% "X1") %>%
  filter(prop >= 0.5 & prop <= 0.5) %>%
  mutate(prop = recode(prop, "0.5" = "Missingness proportion = 0.5")) %>%
  mutate(mech = factor(mech, levels = c("MCAR", "MAR"))) %>%
  ggplot() +
  aes(x = method, y = cov, colour = mech, fill = mech) +
  geom_point(shape = 21) +
  guides(colour = "none") +
  scale_colour_manual(values = c("MCAR" = "#f8766d", "MAR" = "#f8766d")) +
  scale_fill_manual(values = c("MCAR" = "gray74", "MAR" = "gray29")) +
  labs(
    # subtitle = 
    # "Coverage of X1 for different imputation methods
    # under MCAR and MAR missingness mechanisms with a missingness proportion of 50%",
    y = "Coverage",
    x = "Imputation method",
    fill = 
      "Missingness mechanism"
  ) +
  theme_bw() +
  theme(legend.position = "top") +
  facet_nested_wrap(vars(prop, term), scales = "fixed") +
  geom_hline(yintercept = 0.94, linetype = "dashed", color = "darkgray") + 
  geom_hline(yintercept = 0.96, linetype = "dashed", color = "darkgray") + 
  geom_hline(yintercept = 0.95, linetype="dotted", color = "lightgray") +
  theme(text = element_text(size = 11,  family = "mono"))

plot(cov_plot_X1)

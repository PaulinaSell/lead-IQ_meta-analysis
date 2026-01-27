# Sensitivity analyses ----

# 1 Risk of Bias ----
# 1.1 Bayesian MA ----
library(brms)
library(tidyverse)
library(tidybayes)
library(ggridges)
library(glue)
library(bayesplot)
library(HDInterval)
library(bayestestR)
library(posterior)

data_full = read.csv("data/study_data_leadIQloss.csv")
# Subsetting full study base for sensitivity analysis: excluding studies with RoB
data_low_medium <- data_full[
  !data_full$author_year %in% c("Crump 2013", "Earl 2016"),
]
data_low <- data_full[
  !data_full$author_year %in%
    c("Crump 2013", "Earl 2016", "Lucchini 2012", "Lucchini 2019"),
]

# Set priors
priors <- c(
  prior(normal(-1, 2), class = Intercept), # overall effect size µ
  prior(normal(0, 2), class = sd, lb = 0)
) # between-study heterogeneity τ


# store data subsets in list for looping
models_data <- list(
  full = data_full,
  low_medium = data_low_medium,
  low = data_low
)


# create empty list for storing models
m.brm_models <- list(
  main = list(),
  prior = list()
)

# Loop through model fits: full study base -> low & medium RoB studies -> low RoB studies
for (data_subset in names(models_data)) {
  data <- models_data[[data_subset]]

  # main model (with weakly informative priors)
  m.brm <- brm(
    beta_ln | se(se_beta_ln) ~ 1 + (1 | author_year),
    data = data,
    prior = priors,
    save_pars = save_pars(all = TRUE),
    chains = 4,
    iter = 4000
  )

  # Sample from prior only
  fitPrior <- brm(
    beta_ln | se(se_beta_ln) ~ 1 + (1 | author_year),
    data = data,
    prior = priors,
    sample_prior = "only",
    chains = 4,
    iter = 4000
  )

  # store models in list (env)
  m.brm_models$main[[data_subset]] <- m.brm
  m.brm_models$prior[[data_subset]] <- fitPrior

  # save models for later use
  saveRDS(
    m.brm,
    file = paste0("models/sensitivity/RoB/m.brm_", data_subset, ".rds")
  )
  saveRDS(
    fitPrior,
    file = paste0("models/sensitivity/RoB/fitPrior_", data_subset, ".rds")
  )
  saveRDS(m.brm_models, file = "models/sensitivity/RoB/m.brm_models_list.rds") # save whole list

  # save draws from main models (not fitPrior) as .csv for EBD assessment (later)
  draws <- spread_draws(m.brm, b_Intercept, sd_author_year__Intercept)
  write.csv(
    draws,
    paste0("results/draws_beta_tau_", data_subset, ".csv"),
    row.names = F
  )
}

# 1.2 Frequentist MA ----

library(metafor)

# create list to store results
rma_models <- list()

for (data_subset in names(models_data)) {
  data <- models_data[[data_subset]]

  # run model
  rma_model <- rma(yi = beta_ln, sei = se_beta_ln, data = data, method = "REML") # Restricted Maximum Likelihood

  # save in env
  rma_models[[data_subset]] <- rma_model

  # print model summaries
  cat("\n===========================================\n")
  cat("Model including:", data_subset, "RoB studies \n")
  cat("===========================================\n")
  print(summary(rma_model))

  # save models for later use
  saveRDS(
    rma_model,
    file = paste0("models/sensitivity/RoB/freq_", data_subset, ".rds")
  )
}

# 2 Priors (only for BMA)----

# priors for main analysis:
priors <- c(
  prior(normal(-1, 2), class = Intercept), # overall effect size µ
  prior(normal(0, 2), class = sd, lb = 0)
) # between-study heterogeneity τ

# less informative / wider priors:
priors_wide <- c(
  prior(normal(0, 6), class = Intercept), # overall effect size µ
  prior(normal(0, 6), class = sd, lb = 0)
) # between-study heterogeneity τ

# more informative / more narrow priors:
priors_narrow <- c(
  prior(normal(-1, 1), class = Intercept), # overall effect size µ
  prior(normal(0, 1), class = sd, lb = 0)
) # between-study heterogeneity τ

# store priors in list
priors_list <- list(
  main = priors,
  wide = priors_wide,
  narrow = priors_narrow
)

# loop
for (prior_set in names(priors_list)) {
  priors = priors_list[[prior_set]]

  m.brm <- brm(
    beta_ln | se(se_beta_ln) ~ 1 + (1 | author_year),
    data = data_full,
    prior = priors,
    save_pars = save_pars(all = TRUE),
    chains = 4,
    iter = 4000
  )

  # Sample from prior only
  fitPrior <- brm(
    beta_ln | se(se_beta_ln) ~ 1 + (1 | author_year),
    data = data_full,
    prior = priors,
    sample_prior = "only",
    chains = 4,
    iter = 4000
  )

  # store models in list (env)
  m.brm_models$main[[prior_set]] <- m.brm
  m.brm_models$prior[[prior_set]] <- fitPrior

  # save models for later use
  saveRDS(
    m.brm,
    file = paste0("models/sensitivity/priors/m.brm_", prior_set, ".rds")
  )
  saveRDS(
    fitPrior,
    file = paste0("models/sensitivity/priors/fitPrior_", prior_set, ".rds")
  )
  saveRDS(
    m.brm_models,
    file = "models/sensitivity/priors/m.brm_models_list.rds"
  ) # save whole list
}


# 2.1 Checking Results ----
results_priors <- list(
  m.brm_main = m.brm_main,
  m.brm_wide = m.brm_wide,
  m.brm_narrow = m.brm_narrow,
  fitPrior_main = fitPrior_main,
  fitPrior_wide = fitPrior_wide,
  fitPrior_narrow = fitPrior_narrow
)

for (result_name in names(results_priors)) {
  result = results_priors[[result_name]]

  pp_check_priors <- pp_check(result, ndraws = 20) +
    ggtitle(paste("PP Check:", result_name))
  print(pp_check_priors)

  plot_chains <- plot(
    result,
    variable = c("b_Intercept", "sd_author_year__Intercept")
  ) +
    ggtitle(paste("MCMC Chains for ", result_name, "Model"))
  print(plot_chains)

  # print model summaries
  cat(
    "\n======================================================================================\n"
  )
  cat("Sensitivity analysis with alternative priors. Model:", result_name, "\n")
  cat(
    "======================================================================================\n"
  )
  print(summary(result))
}

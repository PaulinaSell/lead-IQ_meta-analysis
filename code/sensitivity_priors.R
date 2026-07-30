# Sensitivity analyses varying one prior at a time.
# Run from the repository root: Rscript code/sensitivity_priors.R

library(brms)
library(dplyr)
library(tibble)
library(purrr)

source("code/prior_configs.R")

data_full <- read.csv("data/study_data_leadIQloss.csv")

dir.create(ONEWAY_MODEL_DIR, recursive = TRUE, showWarnings = FALSE)

MODEL_FORMULA <- bf(beta_ln | se(se_beta_ln) ~ 1 + (1 | author_year))

fit_config <- function(config, config_name) {
  if (!config$fit) {
    message("Reusing existing fit for ", config_name)
    return(invisible(NULL))
  }

  message("Fitting ", config_name)

  m.brm <- brm(
    MODEL_FORMULA,
    data = data_full,
    prior = config$priors,
    chains = 4,
    iter = 4000,
    control = list(adapt_delta = 0.99)
  )

  fitPrior <- brm(
    MODEL_FORMULA,
    data = data_full,
    prior = config$priors,
    sample_prior = "only",
    chains = 4,
    iter = 4000,
    control = list(adapt_delta = 0.99)
  )

  saveRDS(m.brm, config$posterior_path)
  saveRDS(fitPrior, config$prior_path)

  invisible(NULL)
}

iwalk(oneway_prior_configs, fit_config)

extract_pooled <- function(config, config_name) {
  fe <- fixef(load_posterior(config))

  tibble(
    sensitivity = config$arm,
    variant = config_name,
    variant_label = config$label,
    model_type = "Bayesian",
    mean = fe["Intercept", "Estimate"],
    lb = fe["Intercept", "Q2.5"],
    ub = fe["Intercept", "Q97.5"]
  )
}

table_sensi_oneway <- imap_dfr(oneway_prior_configs, extract_pooled)

saveRDS(
  table_sensi_oneway,
  "manuscript/tables_figures/main/table_sensi_oneway.rds"
)

print(as.data.frame(table_sensi_oneway))

# Sensitivity analyses ----

# 1 Risk of Bias ----
# 1.1 Bayesian MA ----
library(brms)
library(metafor)
library(tidyverse)
library(tidybayes)
library(ggridges)
library(glue)
library(bayesplot)
library(HDInterval)
library(bayestestR)
library(posterior)

data_full <- read.csv("data/study_data_leadIQloss.csv")
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

# Loop through brms model fits
# Full study base -> low & medium RoB studies -> low RoB studies
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
  m.brm_models$rob$main[[data_subset]] <- m.brm
  m.brm_models$rob$prior[[data_subset]] <- fitPrior

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
  priors <- priors_list[[prior_set]]

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
  m.brm_models$prior$main[[prior_set]] <- m.brm
  m.brm_models$prior$prior[[prior_set]] <- fitPrior

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
  )
}

# read models
# m.brm_models <- readRDS("models/sensitivity/priors/m.brm_models_list.rds")

# 2.1 Checking Results ----
# results_priors <- list(
#   m.brm_main = m.brm_models$main$main,
#   m.brm_wide = m.brm_models$main$wide,
#   m.brm_narrow = m.brm_models$main$narrow,
#   fitPrior_main = m.brm_models$prior$main,
#   fitPrior_wide = m.brm_models$prior$wide,
#   fitPrior_narrow = m.brm_models$prior$narrow
# )

# for (result_name in names(results_priors)) {
#   result = results_priors[[result_name]]

#   pp_check_priors <- pp_check(result, ndraws = 20) +
#     ggtitle(paste("PP Check:", result_name))
#   print(pp_check_priors)

#   plot_chains <- plot(
#     result,
#     variable = c("b_Intercept", "sd_author_year__Intercept")
#   )

#   # print model summaries
#   cat(
#     "\n======================================================================================\n"
#   )
#   cat("Sensitivity analysis with alternative priors. Model:", result_name, "\n")
#   cat(
#     "======================================================================================\n"
#   )
#   print(summary(result))
# }

# Table for publication
m.brm_models <- readRDS("models/sensitivity/priors/m.brm_models_list.rds") # both priors and RoB sensi analysis results in here
freq_models <- list(
  full = readRDS("models/sensitivity/RoB/freq_full.rds"),
  low_medium = readRDS("models/sensitivity/RoB/freq_low_medium.rds"),
  low = readRDS("models/sensitivity/RoB/freq_low.rds")
)

library(purrr)

extract_brms <- function(fit, coeff = "Intercept") {
  fe <- fixef(fit)

  tibble(
    bayes_mean = fe[coeff, "Estimate"],
    bayes_lb = fe[coeff, "Q2.5"],
    bayes_ub = fe[coeff, "Q97.5"]
  )
}

extract_rma <- function(fit) {
  stopifnot(inherits(fit, "rma"))

  tibble(
    mean = as.numeric(fit$b[1, 1]),
    lb = fit$ci.lb[1],
    ub = fit$ci.ub[1]
  )
}

# Bayes sensitivity analysis table
table_sensi <- imap_dfr(
  m.brm_models,
  function(sensitivity_group, group_name) {
    posterior_models <- sensitivity_group$main # excluding sample_prior = "only" model fits

    imap_dfr(
      posterior_models,
      function(fit, variant) {
        extract_brms(fit) %>%
          mutate(
            sensitivity_group = group_name,
            sensitivity_variant = variant
          )
      }
    )
  }
)

tab_bayes_rob <- imap_dfr(
  m.brm_models$rob$main,
  function(fit, variant) {
    extract_brms(fit) %>%
      mutate(
        sensitivity = "RoB",
        variant = variant,
        model_type = "Bayesian"
      )
  }
)

tab_bayes_prior <- imap_dfr(
  m.brm_models$prior$main,
  ~ extract_brms(.x) %>%
    mutate(
      sensitivity = "Priors",
      variant = .y,
      model_type = "Bayesian"
    )
)

tab_bayes <- bind_rows(tab_bayes_rob, tab_bayes_prior)

# Freq sensitivity analysis table
tab_freq_rob <- map_dfr(
  names(freq_models),
  function(variant) {
    fit <- freq_models[[variant]]

    extract_rma(fit) %>%
      mutate(
        sensitivity = "RoB",
        variant = variant,
        model_type = "Frequentist"
      )
  }
)


table_all <- bind_rows(
  tab_bayes %>%
    transmute(
      sensitivity,
      variant,
      model_type,
      mean = bayes_mean,
      lb = bayes_lb,
      ub = bayes_ub
    ),
  tab_freq_rob
)

saveRDS(
  table_all,
  file = "manuscript/tables_figures/main/table_sensi.rds"
)

table_pub <- table_all %>%
  mutate(
    Sensitivity = factor(sensitivity, levels = c("Priors", "RoB")),
    Variant = variant,
    "Model type" = factor(model_type, levels = c("Bayesian", "Frequentist")),
    Estimate = sprintf("%.2f (%.2f; %.2f)", mean, lb, ub)
  ) %>%
  arrange(Sensitivity, Variant, "Model type") %>%
  select("Model type", Sensitivity, Variant, Estimate)

knitr::kable(
  table_pub,
  caption = "Sensitivity analyses by modeling framework (Bayesian vs Frequentist)"
)

# Leave-One-Out (loo) cross-validation ----

# Frequentist model
freq_model <- readRDS("models/main/freq_full.rds")
loo_freq <- leave1out(freq_model)
print(loo_freq)

# Bayesian model - manual leave1out equivalent
m.brm <- readRDS("models/main/m.brm_full.rds")

priors <- c(
  prior(normal(-1, 2), class = Intercept), # overall effect size µ
  prior(normal(0, 2), class = sd, lb = 0)
)

n_studies <- nrow(data_full)

loo_bayes_results <- data.frame(
  excluded_study = character(n_studies),
  estimate = numeric(n_studies),
  lower_ci = numeric(n_studies),
  upper_ci = numeric(n_studies)
)

for (i in 1:n_studies) {
  # remove study i from the data
  data_without_study <- data_full[-i, ]

  # refit the model on the remaining 11 studies
  fit_without_study <- brm(
    beta_ln | se(se_beta_ln) ~ 1 + (1 | author_year),
    data = data_without_study,
    prior = priors,
    save_pars = save_pars(all = TRUE),
    chains = 4,
    iter = 4000
  )

  # Pull out the pooled estimate and CI for this fit
  pooled_estimate <- fixef(fit_without_study)["Intercept", ]

  # Store results in table, row by row
  loo_bayes_results$excluded_study[i] <- data_full$author_year[i]
  loo_bayes_results$estimate[i] <- pooled_estimate["Estimate"]
  loo_bayes_results$lower_ci[i] <- pooled_estimate["Q2.5"]
  loo_bayes_results$upper_ci[i] <- pooled_estimate["Q97.5"]

  # save fits for inspection
  saveRDS(
    fit_without_study,
    file = paste0("models/sensitivity/loo/m.brm_excl_", i, ".rds")
  )

  # print progress
  cat("Finished fitting excluding:", data_full$author_year[i], "\n")
}

saveRDS(loo_bayes_results, "results/loo_bayes_results.rds")
saveRDS(loo_freq, "results/loo_freq_results.rds")

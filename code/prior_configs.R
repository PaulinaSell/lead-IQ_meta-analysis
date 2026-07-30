# Prior specifications for the one-at-a-time prior sensitivity analyses.
# Sourced by code/sensitivity_priors.R and code/figures_priors_oneway.R.

library(brms)

ONEWAY_MODEL_DIR <- "models/sensitivity/priors_oneway"

MU_DEFAULT_LOCATION <- -1
MU_SCALE <- 2
MU_FLAT_SCALE <- 100
TAU_DEFAULT_SCALE <- 2

mu_prior <- function(location, scale = MU_SCALE) {
  prior_string(
    sprintf("normal(%.6f, %.6f)", location, scale),
    class = "Intercept"
  )
}

tau_prior <- function(scale) {
  prior_string(sprintf("normal(0, %.6f)", scale), class = "sd", lb = 0)
}

oneway_path <- function(prefix, config_name) {
  file.path(ONEWAY_MODEL_DIR, paste0(prefix, "_", config_name, ".rds"))
}

oneway_prior_configs <- list(
  mu_skeptical = list(
    arm = "Effect prior",
    label = "Skeptical: Normal(-0.33, 2)",
    priors = c(
      mu_prior(MU_DEFAULT_LOCATION / 3),
      tau_prior(TAU_DEFAULT_SCALE)
    ),
    posterior_path = oneway_path("m.brm", "mu_skeptical"),
    prior_path = oneway_path("fitPrior", "mu_skeptical"),
    fit = TRUE
  ),
  mu_default = list(
    arm = "Effect prior",
    label = "Default: Normal(-1, 2)",
    priors = c(
      mu_prior(MU_DEFAULT_LOCATION),
      tau_prior(TAU_DEFAULT_SCALE)
    ),
    posterior_path = "models/main/m.brm_full.rds",
    prior_path = "models/main/fitPrior_full.rds",
    fit = FALSE
  ),
  mu_enthusiastic = list(
    arm = "Effect prior",
    label = "Enthusiastic: Normal(-3, 2)",
    priors = c(
      mu_prior(MU_DEFAULT_LOCATION * 3),
      tau_prior(TAU_DEFAULT_SCALE)
    ),
    posterior_path = oneway_path("m.brm", "mu_enthusiastic"),
    prior_path = oneway_path("fitPrior", "mu_enthusiastic"),
    fit = TRUE
  ),
  tau_strong = list(
    arm = "Heterogeneity prior",
    label = "Strong shrinkage: Half-normal(0.67)",
    priors = c(
      mu_prior(0, MU_FLAT_SCALE),
      tau_prior(TAU_DEFAULT_SCALE / 3)
    ),
    posterior_path = oneway_path("m.brm", "tau_strong"),
    prior_path = oneway_path("fitPrior", "tau_strong"),
    fit = TRUE
  ),
  tau_default = list(
    arm = "Heterogeneity prior",
    label = "Default: Half-normal(2)",
    priors = c(
      mu_prior(0, MU_FLAT_SCALE),
      tau_prior(TAU_DEFAULT_SCALE)
    ),
    posterior_path = oneway_path("m.brm", "tau_default"),
    prior_path = oneway_path("fitPrior", "tau_default"),
    fit = TRUE
  ),
  tau_weak = list(
    arm = "Heterogeneity prior",
    label = "Weak shrinkage: Half-normal(6)",
    priors = c(
      mu_prior(0, MU_FLAT_SCALE),
      tau_prior(TAU_DEFAULT_SCALE * 3)
    ),
    posterior_path = oneway_path("m.brm", "tau_weak"),
    prior_path = oneway_path("fitPrior", "tau_weak"),
    fit = TRUE
  )
)

# Existing analyses that move both priors at once, kept for figure 5d.
joint_prior_configs <- list(
  main = list(
    label = "Main: Normal(-1, 2) and Half-normal(2)",
    posterior_path = "models/main/m.brm_full.rds",
    prior_path = "models/main/fitPrior_full.rds"
  ),
  narrow = list(
    label = "Narrow: Normal(-1, 1) and Half-normal(1)",
    posterior_path = "models/sensitivity/priors/m.brm_narrow.rds",
    prior_path = "models/sensitivity/priors/fitPrior_narrow.rds"
  ),
  wide = list(
    label = "Wide: Normal(0, 6) and Half-normal(6)",
    posterior_path = "models/sensitivity/priors/m.brm_wide.rds",
    prior_path = "models/sensitivity/priors/fitPrior_wide.rds"
  )
)

load_posterior <- function(config) readRDS(config$posterior_path)

load_prior <- function(config) readRDS(config$prior_path)

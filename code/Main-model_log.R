# Bayesian Meta Analysis Lead - IQ loss in Children
# Hierarchical Model: Lvl 1: individual participants; Lvl 2: participants nested within studies

library(brms)
library(tidybayes)
library(tidyverse)

data = read.csv("data/study_data_leadIQloss.csv")

# Set priors ----
priors <- c(
  prior(normal(-1, 2), class = Intercept), # overall effect size µ
  prior(normal(0, 2), class = sd, lb = 0)
) # between-study heterogeneity τ


# Main model
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

# save models
# saveRDS(m.brm, "models/main/m.brm_full.rds")
# saveRDS(fitPrior, file = "models/main/fitPrior_full.rds")

# read models
m.brm <- readRDS("models/main/m.brm_full.rds")
fitPrior <- readRDS("models/main/fitPrior_full.rds")

# save draws from main models (not fitPrior) as .csv for EBD assessment (later)
draws <- spread_draws(m.brm, b_Intercept, sd_author_year__Intercept)
write.csv(
  draws,
  "results/draws_beta_tau_full.csv",
  row.names = F
)

# save draws in EBD folder for later use
write.csv(
  draws,
  "data/draws_beta_tau.csv",
  row.names = F
)

# plot the MCMC chains & posterior distributions
dev.new()
plot(m.brm, variable = c("b_Intercept", "sd_author_year__Intercept"))

summary(m.brm)

pp_check(fitPrior, ndraws = 20)
pp_check(m.brm, ndraws = 20)

# investigate model fit
loo(m.brm, moment_match = TRUE, reloo = TRUE)

# Check ESS & MCSE ----
draws <- as_draws(m.brm)
summarise_draws(draws, "ess_bulk", "ess_tail") # bulk ESS ~1,000 and tail ESS ~2,000 indicate well-mixed chains and reasonably precise estimates
summarise_draws(draws, "mcse_mean", "mcse_sd") # relative MCSE: MCSE / posterior SD = 0.0159/1.91 = 0.0083 -> relative MCSE < 0.01 → excellent (very high precision)


# extract draws for EBD assessment, using spread_draws ----
posterior_summary(m.brm)
draws_pooled_b_sd <- spread_draws(
  m.brm,
  b_Intercept,
  sd_author_year__Intercept
)


write.csv(
  draws_pooled_b_sd,
  "results/draws_pooled_b_sd_logBLL.csv",
  row.names = F
)

# calculate I2 from tau
# Extract posterior samples of tau
tau_samples <- as_draws_df(m.brm, variable = "sd_author_year__Intercept")
tau_sq_samples <- tau_samples$sd_author_year__Intercept^2

# Get the SAME typical within-study variance that metafor used
freq_model <- readRDS(
  "models/main/freq_full.rds"
)
vt <- freq_model$vt

# Calculate I2 using metafor's approach for each posterior draw
I2_samples <- tau_sq_samples / (tau_sq_samples + vt)

# Summarize tau2
tau_sq_mean <- mean(tau_sq_samples)
tau_sq_median <- median(tau_sq_samples)
tau_sq_ci <- quantile(tau_sq_samples, c(0.025, 0.975))

# Summarize I2
I2_posterior_mean <- mean(I2_samples)
I2_posterior_median <- median(I2_samples)
I2_credible_interval <- quantile(I2_samples, c(0.025, 0.975))

# summary for Rhat
s_brm <- summary(m.brm)

# Calculate PPI for hypothetical new studies
mu <- draws$b_Intercept
tau <- draws$sd_author_year__Intercept
theta_new <- rnorm(nrow(draws), mu, tau)
theta_new_mean <- mean(theta_new)
theta_new_qs <- quantile(theta_new, probs = c(0.025, 0.975))

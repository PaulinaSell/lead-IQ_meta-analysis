# Bayesian Meta Analysis Lead - IQ loss in Children
# Hierarchical Model: Lvl 1: individual participants; Lvl 2: participants nested within studies

rm(list = ls(all = T))

library(brms)
library(tidyverse)
library(tidybayes)
library(ggridges)
library(glue)
library(bayesplot)
library(HDInterval)
library(bayestestR)
library(posterior)

data = read.csv("data/study_data_leadIQloss.csv")

# Set priors ----
priors <- c(
  prior(normal(-1, 2), class = Intercept), # overall effect size µ
  prior(normal(0, 2), class = sd, lb = 0)
) # between-study heterogeneity τ‚


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

# save models for later use
saveRDS(m.brm, "models/main/m.brm_full.rds")
saveRDS(fitPrior, file = "models/main/fitPrior_full.rds")

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
  "/Users/paulinasell/Documents/UBA/PARC/R/EBD Lead - IQ loss/Project_lead-IQloss/data/draws_beta_tau.csv",
  row.names = F
)

# plot the MCMC chains & posterior distributions
dev.new()
plot(m.brm, variable = c("b_Intercept", "sd_author_year__Intercept"))

summary(m.brm)

pp_check(fitPrior, ndraws = 20)
pp_check(m.brm, ndraws = 20)

# investigate model fit
loo(m.brm, moment_match = TRUE, reloo = TRUE) # Leave-One-Out Cross-Validation (LOO-CV)
# reloo = T because 1 approximation was still bad, suggested by output:
# We recommend to set 'reloo = TRUE' in order to calculate the ELPD without the assumption that these observations are negligible.
# This will refit the model 1 times to compute the ELPDs for the problematic observations directly.

# Check ESS & MCSE ----

# ESS (bulk) = effective sample size for the bulk of the posterior (i.e. central mass).
# It answers: “How many independent draws are these autocorrelated MCMC samples worth for estimating central summaries like the mean/median?”
# ESS (tail) = effective sample size for the tails of the posterior (useful for accurate extreme quantiles, credible intervals, tail probabilities).

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

# change names before exporting again!
# write.csv(draws_pooled_b_sd, "/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/results/draws_pooled_b_sd_logBLL_12studies.csv", row.names = F)
# write.csv(draws_pooled_b_sd, "/Users/paulinasell/Documents/UBA/PARC/R/EBD Lead - IQ loss/Project_lead-IQloss/data/draws_pooled_b_sd_logBLL_12studies.csv", row.names = F)

# calculate I2 from tau

# 1. Extract tau (between-study SD)
posterior_summary(m.brm, variable = "sd_author_year__Intercept")
# This gives you the posterior mean, SD, and credible intervals

# 2. Get tau-squared from tau
tau_samples <- as_draws_df(m.brm, variable = "sd_author_year__Intercept")
tau_sq_samples <- tau_samples$sd_author_year__Intercept^2

# Posterior summary of tau-squared
mean(tau_sq_samples) # Compare to 3.7946
quantile(tau_sq_samples, c(0.025, 0.975)) # 95% CrI

# 3. Calculate I² for each study and each posterior draw
n_studies <- nrow(data)
n_samples <- length(tau_sq_samples)

I2_by_study <- matrix(NA, nrow = n_samples, ncol = n_studies)

data$vi <- data$se_beta_ln^2

for (i in 1:n_studies) {
  I2_by_study[, i] <- tau_sq_samples / (tau_sq_samples + data$vi[i])
}

# 4. Average I² across studies for each posterior draw
I2_samples <- rowMeans(I2_by_study)

# 5. Summarize the posterior distribution
I2_posterior_mean <- mean(I2_samples)
I2_posterior_median <- median(I2_samples)
I2_credible_interval <- quantile(I2_samples, c(0.025, 0.975))

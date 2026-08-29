# Frequentist Random Effects Meta Analysis Lead - IQ loss in Children

library(metafor)
library(tidyverse)

data_full = read.csv("data/study_data_leadIQloss.csv")

# run unadjusted model
rma_model <- rma(
  yi = beta_ln,
  sei = se_beta_ln,
  data = data_full,
  method = "REML"
) # Restricted Maximum Likelihood

print(rma_model)

rma_model_knha <- rma(
  yi = beta_ln,
  sei = se_beta_ln,
  data = data_full,
  method = "REML",
  test = "knha"
) # Small sample adjustment

print(rma_model_knha)

saveRDS(rma_model_knha, file = "models/main/freq_full.rds")

# Funnel plot
funnel(rma_model_knha)

# Diagnostics
# Influence diagnostics (equivalent to Pareto-k checks)
inf <- influence(rma_model_knha)
plot(inf)

# Residuals
plot(rma_model_knha)

# Leave-One-Out Analysis
# loo_results <- leave1out(rma_model_knha)
# print(loo_results)

forest(rma_model_knha)

# Parametric simulation of predictive intervals (hypothetical new studies)
# 8000 draws to match the Bayesian analysis
n_samples <- 8000
mu_hat <- as.numeric(coef(rma_model_knha))
se_mu <- rma_model_knha$se
tau_hat <- sqrt(rma_model_knha$tau2)
se_tau <- rma_model_knha$se.tau2 / (2 * tau_hat)

mu_samples <- rnorm(n_samples, mu_hat, se_mu)
log_tau_mu <- log(tau_hat)
log_tau_sd <- se_tau / tau_hat
tau_samples <- exp(rnorm(n_samples, log_tau_mu, log_tau_sd))
theta_new <- rnorm(n_samples, mu_samples, tau_samples)
theta_new_mean <- mean(theta_new)
theta_new_qs <- quantile(theta_new, probs = c(0.025, 0.975))

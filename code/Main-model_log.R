# Bayesian Meta Analysis Lead - IQ loss in Children
# Hierarchical Model: Lvl 1: individual participants; Lvl 2: participants nested within studies

rm(list=ls(all=T))

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
priors <- c(prior(normal(-1, 2), class = Intercept), # overall effect size µ
             prior(normal(0, 2), class = sd, lb = 0)) # between-study heterogeneity τ‚


# main model (with weakly informative priors) and full dataset
m.brm <- brm(
    beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
    data = data,
    prior = priors,
    save_pars = save_pars(all = TRUE),
    chains = 4,
    iter = 4000)
  
  # Sample from prior only
fitPrior <- brm(
    beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), 
    data = data, 
    prior = priors,
    sample_prior = "only",
    chains = 4,
    iter = 4000)

# save models for later use
saveRDS(m.brm, "models/main/m.brm_full.rds")
saveRDS(fitPrior, file = "models/main/fitPrior_full.rds")

# save draws from main models (not fitPrior) as .csv for EBD assessment (later)
draws <- spread_draws(m.brm, b_Intercept, sd_author_year__Intercept)
write.csv(draws, paste0("results/draws_beta_tau_", data_subset, ".csv"), row.names = F)


# plot the MCMC chains & posterior distributions
plot(m.brm_full, variable = c("b_Intercept", "sd_author_year__Intercept"))

summary(m.brm_full)

# Sample from prior only
fitPrior_full <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), 
  data = data_full, 
  prior = priors,
  sample_prior = "only",
  chains = 4,
  iter = 4000)

# save model for manuscript
# saveRDS(fitPrior_full, file = "models/fitPrior_main_full.rds")

plot(fitPrior_full, variable = c("b_Intercept", "sd_author_year__Intercept"))

pp_check(fitPrior_full, ndraws = 20)
pp_check(m.brm_full, ndraws = 20)

# investigate model fit
loo(m.brm_full,
    moment_match = TRUE, 
    reloo = TRUE) # Leave-One-Out Cross-Validation (LOO-CV)
                                              # reloo = T because 1 approximation was still bad, suggested by output: 
                                                    # We recommend to set 'reloo = TRUE' in order to calculate the ELPD without the assumption that these observations are negligible. 
                                                    # This will refit the model 1 times to compute the ELPDs for the problematic observations directly. 




# Check ESS & MCSE ----

# ESS (bulk) = effective sample size for the bulk of the posterior (i.e. central mass). 
# It answers: “How many independent draws are these autocorrelated MCMC samples worth for estimating central summaries like the mean/median?”
# ESS (tail) = effective sample size for the tails of the posterior (useful for accurate extreme quantiles, credible intervals, tail probabilities).

draws_full <- as_draws(m.brm_full)
summarise_draws(draws_full, "ess_bulk", "ess_tail") # bulk ESS ~1,000 and tail ESS ~2,000 indicate well-mixed chains and reasonably precise estimates
summarise_draws(draws_full, "mcse_mean", "mcse_sd") # relative MCSE: MCSE / posterior SD = 0.0159/1.91 = 0.0083 -> relative MCSE < 0.01 → excellent (very high precision)


# extract draws for EBD assessment, using spread_draws ----
posterior_summary(m.brm_full)
draws_pooled_b_sd <- spread_draws(m.brm_full, b_Intercept, sd_author_year__Intercept)

# change names before exporting again!
# write.csv(draws_pooled_b_sd, "/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/results/draws_pooled_b_sd_logBLL_12studies.csv", row.names = F)
# write.csv(draws_pooled_b_sd, "/Users/paulinasell/Documents/UBA/PARC/R/EBD Lead - IQ loss/Project_lead-IQloss/data/draws_pooled_b_sd_logBLL_12studies.csv", row.names = F)

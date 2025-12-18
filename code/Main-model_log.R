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

data_full = read.csv("/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/data/study_data_leadIQloss.csv")

# Subsetting full study base for sensitivity analysis: excluding studies with RoB
data_low_medium <- data_full[!data_full$author_year %in% c("Crump 2013", "Earl 2016"), ]
data_low <- data_full[!data_full$author_year %in% c("Crump 2013", "Earl 2016", "Lucchini 2012", "Lucchini 2019"), ]

# Set priors ----
priors <- c(prior(normal(-1, 2), class = Intercept), # overall effect size µ
             prior(normal(0, 2), class = sd, lb = 0)) # between-study heterogeneity τ


# Fit model with full study base (12 studies) ----

# Main model (with weakly informative priors)
m.brm_full <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data_full,
  prior = priors,
  # control = list(adapt_delta = 0.90),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000)

# save model
saveRDS(m.brm_full, file = "models/m.brm_main_full.rds")

# plot the MCMC chains & posterior distributions
plot(m.brm_full, variable = c("b_Intercept", "sd_author_year__Intercept"))

summary(m.brm_full)

# Sample from prior only
fitPrior_full <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), 
  data = data_full, 
  prior = priors,
  # control = list(adapt_delta = 0.90),
  sample_prior = "only",
  chains = 4,
  iter = 4000)

# save model for manuscript
saveRDS(fitPrior_full, file = "models/fitPrior_main_full.rds")

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

# Sensitivity analysis 1: excluding studies with high RoB (n=2), leaving 10 studies ----

# Main model (with weakly informative priors)
m.brm_low_medium <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data_low_medium,
  prior = priors,
  # control = list(adapt_delta = 0.90),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000)

# save model
saveRDS(m.brm_low_medium, file = "models/m.brm_main_low_medium.rds")

# plot the MCMC chains & posterior distributions
plot(m.brm_low_medium, variable = c("b_Intercept", "sd_author_year__Intercept"))

summary(m.brm_low_medium)

# Sample from prior only
fitPrior_low_medium <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), 
  data = data_low_medium, 
  prior = priors,
  # control = list(adapt_delta = 0.90),
  sample_prior = "only",
  chains = 4,
  iter = 4000)

# save model for manuscript
saveRDS(fitPrior_low_medium, file = "models/fitPrior_main_low_medium.rds")

plot(fitPrior_low_medium, variable = c("b_Intercept", "sd_author_year__Intercept"))

pp_check(fitPrior_low_medium, ndraws = 20)
pp_check(m.brm_low_medium, ndraws = 20)

# investigate model fit
loo(m.brm_low_medium,
    moment_match = TRUE, 
    reloo = TRUE)


# Sensitivity analysis 2: excluding studies with high OR moderate RoB (n=4), leaving 8 studies ----

# Main model (with weakly informative priors)
m.brm_low <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data_low,
  prior = priors,
  control = list(adapt_delta = 0.90),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000)

# save model
saveRDS(m.brm_low, file = "models/m.brm_main_low.rds")

# plot the MCMC chains & posterior distributions
plot(m.brm_low, variable = c("b_Intercept", "sd_author_year__Intercept"))

summary(m.brm_low)

# Sample from prior only
fitPrior_low <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), 
  data = data_low, 
  prior = priors,
  control = list(adapt_delta = 0.90),
  sample_prior = "only",
  chains = 4,
  iter = 4000)

# save model for manuscript
saveRDS(fitPrior_low, file = "models/fitPrior_main_low.rds")

plot(fitPrior_low, variable = c("b_Intercept", "sd_author_year__Intercept"))

pp_check(fitPrior_low, ndraws = 20)
pp_check(m.brm_low, ndraws = 20)

# investigate model fit
loo(m.brm_low,
    moment_match = TRUE, 
    reloo = TRUE)



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

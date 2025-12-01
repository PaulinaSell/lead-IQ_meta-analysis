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

data = read.csv("/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/data/study_data_leadIQloss.csv")

# data cleaning: exclude studies that were transformed from linear to log since transformation may not be valid here
data <- data[!data$author_year %in% c("Halabicky 2022", "Iglesias 2011", "Min 2009"), ]


# Set priors ----
priors <- c(prior(normal(-1, 2), class = Intercept), # overall effect size µ
             prior(normal(0, 2), class = sd, lb = 0)) # between-study heterogeneity τ

# Fit model ----

# Main model (with weakly informative priors)
m.brm <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = priors,
  # control = list(adapt_delta = 0.90),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000)

# plot the MCMC chains & posterior distributions
plot(m.brm, variable = c("b_Intercept", "sd_author_year__Intercept"))

summary(m.brm)

# how to extract only mean beta and its sd? 
Bayes_b_Intercept_mean_1 <- round(posterior_summary(m.brm, variable = "b_Intercept")[1], 2)
Bayes_b_Intercept_se_1 <- round(posterior_summary(m.brm, variable = "b_Intercept")[2], 2)

# Sample from prior only
fitPrior <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), 
  data = data, 
  prior = priors,
  # control = list(adapt_delta = 0.90),
  sample_prior = "only",
  chains = 4,
  iter = 4000)

plot(fitPrior, variable = c("b_Intercept", "sd_author_year__Intercept"))

pp_check(fitPrior, ndraws = 20)
pp_check(m.brm, ndraws = 20)

# investigate model fit
loo(m.brm,
    moment_match = TRUE, 
    reloo = TRUE
    ) # Leave-One-Out Cross-Validation (LOO-CV)
                                              # reloo = T because 1 approximation was still bad, suggested by output: 
                                                    # We recommend to set 'reloo = TRUE' in order to calculate the ELPD without the assumption that these observations are negligible. 
                                                    # This will refit the model 1 times to compute the ELPDs for the problematic observations directly. 

# Check Rhat
summary(m.brm)

# Check ESS & MCSE

# ESS (bulk) = effective sample size for the bulk of the posterior (i.e. central mass). 
# It answers: “How many independent draws are these autocorrelated MCMC samples worth for estimating central summaries like the mean/median?”
# ESS (tail) = effective sample size for the tails of the posterior (useful for accurate extreme quantiles, credible intervals, tail probabilities).

draws <- as_draws(m.brm)
summarise_draws(draws, "ess_bulk", "ess_tail") # bulk ESS ~1,000 and tail ESS ~2,000 indicate well-mixed chains and reasonably precise estimates
summarise_draws(draws, "mcse_mean", "mcse_sd") # relative MCSE: MCSE / posterior SD = 0.0159/1.91 = 0.0083 -> relative MCSE < 0.01 → excellent (very high precision)


# exact probability of effect being smaller (in this case: greater) than certain value (-0.45 here) (using empirical cumulative distribution function)
# b.ecdf <- ecdf(post.samples$b_Intercept)
# (1 - b.ecdf(0)) * 100 # -0.45 (95% CI -0.66, -0.24) is the result of Philippes meta analysis

# save model
# saveRDS(m.brm, "models/m.brm") # 31.10.2025

# extract draws for EBD assessment, using spread_draws ----
posterior_summary(m.brm)
draws_pooled_b_sd <- spread_draws(m.brm, b_Intercept, sd_author_year__Intercept)

# write.csv(draws_pooled_b_sd, "/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/results/draws_pooled_b_sd_logBLL_12studies.csv", row.names = F)
# write.csv(draws_pooled_b_sd, "/Users/paulinasell/Documents/UBA/PARC/R/EBD Lead - IQ loss/Project_lead-IQloss/data/draws_pooled_b_sd_logBLL_12studies.csv", row.names = F)

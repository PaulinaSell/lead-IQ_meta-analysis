# Bayesian Meta Analysis Lead - IQ loss in Children
# Hierarchical Model: Lvl 1: individual participants; Lvl 2: participants nested within studies

rm(list=ls(all=T))

library(brms)
library(tidyverse)
library(tidybayes)
library(ggridges)
library(glue)
library(extraDistr)
library(bayesplot)
library(HDInterval)

data = read.csv("/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/data/study_data_leadIQloss.csv")

# excluding studies that were transformed from linear to log, since transformation may not be valid here
data <- data[!data$author_year %in% c("Halabicky 2022", "Iglesias 2011", "Min 2009", "Roy 2013"), ]

# lets look at the data 
ggplot(data, aes(x = 1:nrow(data), y = beta_ln)) +
  geom_hline(yintercept = mean(data$beta_ln), colour = "red") +
  geom_hline(yintercept = median(data$beta_ln), colour = "blue") +
  geom_point() +
  labs(x = "Study", y = "Main Effect Beta") +
  theme_minimal()

ggplot(data, aes(x = 1:nrow(data), y = se_beta_ln)) +
  geom_hline(yintercept = mean(data$se_beta_ln), colour = "red") +
  geom_hline(yintercept = median(data$se_beta_ln), colour = "blue") +
  geom_point() +
  labs(x = "Study", y = "Heterogeneity Tau") +
  theme_minimal()

# how to set priors?
# main effect, beta
# visualize distr
mean <- -2
sd <- 1.5

x <- seq(-8, 4, by = 0.1)
y <- dnorm(x, mean = mean, sd = sd)
ggplot() +
  aes(x, y) +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_line(colour = "blue", linewidth = 2) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ggsave("/Users/paulinasell/Documents/UBA/PARC/Meetings/Meeting with Philippe 04.09.2025/Vis/Prior_mu_norm.png", width = 20, height = 15, units = "cm")

1 - pnorm(0, mean, sd = sd)

# heterogeneity, tau with cauchy
sigma = 0.2
phcauchy(0.25, sigma = sigma) # check probability of less than 0.2 between-study heterogeneity τ in half-cauchy distribution with sigma 0.25 
inverseCDF(c(0.025, 0.975), phcauchy, sigma = sigma) # use inverse Cumulative Density Function to find Q2.5 and Q97.5 of the half-cauchy with sigma x

x <- seq(0, 5, by = 0.01)
y <- dhcauchy(x, sigma = sigma)

ggplot() +
  aes(x, y) +
  geom_line(colour = "orange")

# heterogeneity, tau with inv gamma
shape <- 1.2
scale <- 1

X <- seq(-1, 8, by = 0.1)
Y <- dgamma(X, shape = shape, scale = scale)

ggplot() +
  aes(X, Y) +
  geom_line(colour = "orange")

# heterogeneity, tau with trunc normal
x <- seq(0, 6, by = 0.1)
y <- dnorm(x, mean = 0, sd = 1)
ggplot() +
  aes(x, y) +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_line(colour = "orange", linewidth = 2) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ggsave("/Users/paulinasell/Documents/UBA/PARC/Meetings/Meeting with Philippe 04.09.2025/Vis/Prior_tau_trunc-norm.png", width = 20, height = 15, units = "cm")

# Set priors ----
priors <- c(prior(normal(-1, 2), class = Intercept), # overall effect size µ
             prior(normal(0, 2), class = sd, lb = 0)) # between-study heterogeneity τ

# Fit model ----

# Main model (with weakly informative priors)
m.brm <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = priors,
  # control = list(adapt_delta = 0.99),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000,
  seed = 1223)

# plot the MCMC chains & posterior distributions
plot(m.brm, variable = c("b_Intercept", "sd_author_year__Intercept"))

# # nicer traceplot (incl. warmup)
# posterior_samples_warm = as_draws_df(m.brm, inc_warmup = T)
# names(posterior_samples_warm)[names(posterior_samples_warm) == "b_Intercept"] = "beta"
# names(posterior_samples_warm)[names(posterior_samples_warm) == "sd_author_year__Intercept"] = "sd"

# mcmc_trace(posterior_samples_warm,
#            pars = c("beta", "sd"),
#            facet_args = list(ncol = 1)) +
#   vline_at(2000, color = "red", linetype = 2, size = 0.5)  # mark end of warmup

# Prior predictive check
fitPrior <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), 
  data = data, 
  prior = priors, #1 divergent transitions after warmup
  # control = list(adapt_delta = 0.99),
  sample_prior = "only",
  chains = 4,
  iter = 4000,
  seed = 1223)

pp_check(fitPrior, ndraws = 20)

# Posterior predictive check
pp_check(m.brm, ndraws = 20)

# investigate model fit
loo(m.brm, moment_match = TRUE, reloo = TRUE) # Leave-One-Out Cross-Validation (LOO-CV)
                                              # reloo = T because 1 approximation was still bad, suggested by output: 
                                                    # We recommend to set 'reloo = TRUE' in order to calculate the ELPD without the assumption that these observations are negligible. 
                                                    # This will refit the model 1 times to compute the ELPDs for the problematic observations directly. 


# Check Rhat values & interpret results
summary(m.brm)
ranef(m.brm)


# Trying different priors ----

m.brm1 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(normal(-2, 1), class = Intercept), # overall effect size mu
            prior(normal(0, 1), class = sd, lb = 0)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000,
  seed = 1223)

m.brm2 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(normal(-1, 1), class = Intercept), # overall effect size mu
            prior(normal(0, 1), class = sd, lb = 0)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000,
  seed = 1223)

m.brm3 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(normal(-1, 2), class = Intercept), # overall effect size mu
            prior(normal(0, 1), class = sd, lb = 0)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000,
  seed = 1223)

m.brm4 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(normal(-1, 2), class = Intercept), # overall effect size mu
            prior(normal(0, 2), class = sd, lb = 0)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000,
  seed = 1223)

m.brm5 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(uniform(-5, 2), class = Intercept), # overall effect size mu
            prior(uniform(0, 3), class = sd, lb = 0)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000,
  seed = 1223)

m.brm6 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(student_t(3, -2.9, 2.8), class = Intercept), # overall effect size mu
            prior(student_t(3, 0, 2.8), class = sd)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000,
  seed = 1223)


fitPrior4 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), 
  data = data, 
  prior = c(prior(normal(-1, 2), class = Intercept), # overall effect size mu
            prior(normal(0, 2), class = sd, lb = 0)), #heterogeneity tau
  sample_prior = "only",
  iter = 4000,
  seed = 1223)


summary(m.brm1) # mu normal(-2,1), tau truncnorm(0,1),                  beta mean -1.77 (-2.84    -0.80)
summary(m.brm2) # mu normal(-1,1), tau truncnorm(0,1),                  beta mean -1.49 (-2.50    -0.52)
summary(m.brm3) # mu normal(-1,2), tau truncnorm(0,1),                  beta mean -1.65 (-2.85    -0.58)
summary(m.brm4) # mu normal(-1,2), tau truncnorm(0,2),                  beta mean -1.69 (-3.02    -0.46)
summary(m.brm5) # mu uniform(-5,2), tau uniform(0,3),                   beta mean -1.84 (-3.42    -0.51) 1757 divergent transitions
summary(m.brm6) # mu student_t(3, -2.9, 2.8), tau student_t(3, 0, 2.8), beta mean -1.92 (-3.45    -0.64) 2 divergent transitions        from get_prior

pp_check(fitPrior4, ndraws = 20)
pp_check(m.brm4, ndraws = 20)
loo(m.brm4, moment_match = T)

# saveRDS(m.brm1, "models/m.brm1")
# saveRDS(m.brm2, "models/m.brm2")
# saveRDS(m.brm3, "models/m.brm3")
# saveRDS(m.brm4, "models/m.brm4")
# saveRDS(m.brm5, "models/m.brm5")
# saveRDS(m.brm6, "models/m.brm6")



get_prior(beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), data = data)


# visualize students t distribution
df <- 3
lacation <- 0
scale <- 2.8
x = seq(-5, 5, 0.01)

curve(dt((x - lacation)/scale, df)/scale, from = -5, to = 5, col = "red")



# lets look at posterior distributions ----
post.samples <- as_draws_df(m.brm, variable = c("b_Intercept", "sd_author_year__Intercept"))

# generate density plot for posterior distributions
ggplot(aes(x = b_Intercept), data = post.samples) +
  geom_density(fill = "steelblue",                # set the color
               color = "steelblue", alpha = 0.7) +  
  # geom_vline(xintercept = mean(post.samples$b_Intercept), 
                               # linetype = "dotted", 
                               # color = "red") +
  # geom_vline(xintercept = Mode(post.samples$b_Intercept), 
             # linetype = "dotted", 
             # color = "blue") +
  geom_vline(xintercept = 0,
             color = "grey") +
  labs(x = expression(italic(b_Intercept)),
       y = element_blank()) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# ggsave("/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/results/posterior_dist_b_log_no-Halabicky-Iglesias-Min-Roy_minimal.png", width = 25, height = 15, units = "cm")

ggplot(aes(x = sd_author_year__Intercept), data = post.samples) +
  geom_density(fill = "lightgreen",               # set the color
               color = "lightgreen", alpha = 0.7) +  
  # geom_vline(xintercept = mean(post.samples$sd_author_year__Intercept), 
  #            linetype = "dotted", 
  #            color = "red") +
  # geom_vline(xintercept = Mode(post.samples$sd_author_year__Intercept), 
  #            linetype = "dotted", 
  #            color = "blue") +
  geom_vline(xintercept = 0,
             color = "grey") +
  labs(x = expression(sd_author_year__Intercept),
       y = element_blank()) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

ggsave("/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/results/posterior_dist_sd_log_no-Halabicky-Iglesias-Min-Roy.png", width = 25, height = 15, units = "cm")

# check exact probability of effect being smaller (in this case: greater) than certain value (-0.45 here) (using empirical cumulative distribution function)
b.ecdf <- ecdf(post.samples$b_Intercept)
(1 - b.ecdf(0)) * 100 # -0.45 (95% CI -0.66, -0.24) is the result of Philippes meta analysis



# lets generate a forest plot ----
# first, we have to prepare the data and extract the posterior distribution for each individual study 
study.draws <- spread_draws(m.brm, r_author_year[author_year,], b_Intercept) %>% 
  mutate(b_Intercept = r_author_year + b_Intercept)

# Next: generate distribution of the pooled effect
pooled.effect.draws <- spread_draws(m.brm, b_Intercept) %>% 
  mutate(author_year = "Pooled Effect")

# Next: bind study.draws and pooled.effect.draws to one data frame. 
forest.data <- bind_rows(study.draws, 
                         pooled.effect.draws) %>% 
  ungroup() %>%
  mutate(author_year = str_replace_all(author_year, "[.]", " ")) %>%
  mutate(author_year = fct_rev(factor(author_year)))

# Lastly: forest plot should also display the effect size (SMD and credible interval) of each study. 
# We use forest.data data set, group it by Author, and then use the mean_qi function to calculate these values
forest.data.summary <- group_by(forest.data, author_year) %>% 
  mean_qi(b_Intercept)

ggplot(aes(b_Intercept, 
           relevel(author_year, "Pooled Effect", # Y-axis uses author_year but reorders it so "Pooled Effect" appears at the bottom (after = Inf)
                   after = Inf)), 
       data = forest.data) +
  
  # Box to highlight pooled effect
  geom_rect(data = filter(forest.data.summary, author_year == "Pooled Effect"),
            aes(xmin = -Inf, xmax = Inf, 
                ymin = as.numeric(factor(author_year)) - 0.4, 
                ymax = as.numeric(factor(author_year)) + 0.4),
            fill = "steelblue3", alpha = 0.2) +
  
  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = fixef(m.brm)[1, 1], # line at the pooled effect estimate
             color = "gray50", linewidth = 0.8, linetype = 2) + 
  #geom_vline(xintercept = fixef(m.brm)[1, 3:4], # lines for Q2.5 & Q97.5 from pooled effect (95 % CI)
             #color = "gray50", linewidth = 0.8, linetype = 2) + 
  geom_vline(xintercept = 0, color = "gray20", 
             linewidth = 1) + # line at zero (null effect line)
  
  # Add density ridges
  geom_density_ridges(fill = "steelblue3", # ridge density plots for each study/row
                      rel_min_height = 0.01, # removes very small density values
                      col = NA, scale = 1, # no outline color, scale of distributions
                      alpha = 0.8) + # opacity
  
  # Add point estimates with confidence intervals
  geom_pointinterval(data = forest.data.summary,
                     aes(xmin = .lower, xmax = .upper),
                     fatten_point = 1.5,
                     linewidth = 0.8) +
  
  # Add text and labels
  geom_text(data = mutate_if(forest.data.summary,
                             is.numeric, round, 2), # decimals
            aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"), # format text -> ß [CIlow, CIhigh]
                x = Inf), hjust = "inward") + #aligns text inward from the edge
  labs(x = "beta",
       y = element_blank()) +
  #xlim(-4, 1) +
  theme_light() +
  theme(panel.border = element_blank())

# ggsave("/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/results/forestplot_logBLL_no-Halabicky-Iglesias-Min-Roy.png", width = 25, height = 15, units = "cm")
# ggsave("/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/results/Forest_Priors/m.brm6.png", width = 30, height = 20, units = "cm")

# extract draws for EBD assessment, using spread_draws ----
posterior_summary(m.brm)

draws_pooled_b_sd <- spread_draws(m.brm, b_Intercept, sd_author_year__Intercept)


# write.csv(draws_pooled_b_sd, "/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/results/draws_pooled_b_sd_logBLL_no-Halabicky-Iglesias-Min-Roy.csv", row.names = F)
# write.csv(draws_pooled_b_sd, "/Users/paulinasell/Documents/UBA/PARC/R/EBD Lead - IQ loss/Project_lead-IQloss/data/draws_pooled_b_sd_logBLL_no-Halabicky-Iglesias-Min-Roy.csv", row.names = F)



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

data = read.csv("data/study_data_leadIQloss.csv")

# try excluding studies that were transformed from linear to log, since transformation may not be valid here
data <- data[!data$author_year %in% c("Halabicky 2022", "Iglesias 2011", "Min 2009"), ]

# lets look at the data 
ggplot(data, aes(x = 1:nrow(data), y = beta_ln)) +
  geom_point() +
  labs(x = "Study", y = "Main Effect Beta") +
  theme_minimal()

ggplot(data, aes(x = 1:nrow(data), y = se_beta_ln)) +
  geom_point() +
  labs(x = "Study", y = "Heterogeneity Tau") +
  theme_minimal()

# how to set priors?
# main effect, beta
# visualize distr
x <- seq(-3, 3, by = 0.1)
y <- dnorm(x, mean = 0, sd = 1)
ggplot() +
  aes(x, y) +
  geom_line(colour = "blue")

# heterogeneity, tau
sigma = 0.2
phcauchy(0.25, sigma = sigma) # check probability of less than 0.2 between-study heterogeneity τ in half-cauchy distribution with sigma 0.25 
inverseCDF(c(0.025, 0.975), phcauchy, sigma = sigma) # use inverse Cumulative Density Function to find Q2.5 and Q97.5 of the half-cauchy with sigma x
# visualize distr
x <- seq(0, 5, by = 0.01)
y <- dhcauchy(x, sigma = sigma)
ggplot() +
  aes(x, y) +
  geom_line(colour = "orange")

# Set priors ----
priors <- c(prior(normal(-0.4,0.8), class = Intercept), # overall effect size µ
            prior(cauchy(0,0.2), class = sd, lb = 0)) # between-study heterogeneity τ with lb (lower bound) 0 to make it half-Cauchy

priors_cons <- c(prior(normal(0,1), class = Intercept), 
            prior(cauchy(0,0.5), class = sd, lb = 0)) 


# Fit model ----

# Main model (with conservative priors)
m.brm <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = priors_cons,
  control =                      # Whenever you see the warning "There were x divergent transitions after warmup." you should really think about increasing adapt_delta.
    list(adapt_delta = 0.99),    # Increasing adapt_delta will slow down the sampler but will decrease the number of divergent transitions threatening the validity of your posterior draws
  save_pars = save_pars(all = TRUE),
  iter = 4000,
  seed = 13513)

# plot the MCMC chains & posterior distributions
plot(m.brm) #kommentar


# Prior predictive check
fitPrior <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), 
  data = data, 
  prior = priors_cons,
  sample_prior = "only",
  iter = 4000,
  seed = 13513)

# Plot
pp_check(fitPrior, ndraws = 10)

# investigate model fit
loo(m.brm, moment_match = TRUE, reloo = TRUE) # Leave-One-Out Cross-Validation (LOO-CV)
                                              # reloo = T because 1 approximation was still bad, suggested by output: 
                                                    # We recommend to set 'reloo = TRUE' in order to calculate the ELPD without the assumption that these observations are negligible. 
                                                    # This will refit the model 1 times to compute the ELPDs for the problematic observations directly. 

# Posterior predictive check
pp_check(m.brm, ndraws = 15)

# Examine the pairs() plot to diagnose sampling problems
pairs(m.brm) # especially useful in identifying collinearity between variables (which manifests as narrow bivariate plots) as well as the presence of multiplicative non-identifiabilities (banana-like shapes)


# Check Rhat values & interpret results
summary(m.brm)
ranef(m.brm)

# lets look at posterior distributions ----
# prepare data
str(as_draws_df(m.brm))
post.samples <- as_draws_df(m.brm, variable = c("b_Intercept", "sd_author_year__Intercept"))

# generate density plot for posterior distributions
ggplot(aes(x = b_Intercept), data = post.samples) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_vline(xintercept = mean(post.samples$b_Intercept), 
                               linetype = "dotted", 
                               color = "red") +
  geom_vline(xintercept = Mode(post.samples$b_Intercept), 
             linetype = "dotted", 
             color = "blue") +
  labs(x = expression(italic(b_Intercept)),
       y = element_blank()) +
  theme_minimal()

#ggsave("results/posterior_dist_b_log_no-Halabicky-Iglesias-Min.png", width = 25, height = 15, units = "cm")


ggplot(aes(x = sd_author_year__Intercept), data = post.samples) +
  geom_density(fill = "lightgreen",               # set the color
               color = "lightgreen", alpha = 0.7) +  
  geom_vline(xintercept = mean(post.samples$sd_author_year__Intercept), 
             linetype = "dotted", 
             color = "red") +
  geom_vline(xintercept = Mode(post.samples$sd_author_year__Intercept), 
             linetype = "dotted", 
             color = "blue") +
  labs(x = expression(sd_author_year__Intercept),
       y = element_blank()) +
  theme_minimal()

# ggsave("results/posterior_dist_sd_log_no-Halabicky-Iglesias-Min.png", width = 25, height = 15, units = "cm")


# check exact probability of effect being smaller (in this case: greater) than certain value (-0.45 here) (using empirical cumulative distribution function)
b.ecdf <- ecdf(post.samples$b_Intercept)
b.ecdf(-0.45) # -0.45 (95% CI -0.66, -0.24) is the result of Philippes meta analysis



# lets generate a forest plot ----
# first, we have to prepare the data and extract the posterior distribution for each individual study 
study.draws <- spread_draws(m.brm, r_author_year[author_year,], b_Intercept) %>% 
  mutate(b_Intercept = r_author_year + b_Intercept)

# Next: generate distribution of the pooled effect
pooled.effect.draws <- spread_draws(m.brm, b_Intercept) %>% 
  mutate(author_year = "Pooled Effect")

# Next: bind study.draws and pooled.effect.draws to one data frame. 
# We then start a pipe again, calling ungroup first, and then use mutate to clean the study labels (i.e. replace dots with spaces)
forest.data <- bind_rows(study.draws, 
                         pooled.effect.draws) %>% 
  ungroup() %>%
  mutate(author_year = str_replace_all(author_year, "[.]", " ")) %>%
  mutate(author_year = fct_rev(factor(author_year)))

# Lastly: forest plot should also display the effect size (SMD and credible interval) of each study. 
# We use forest.data data set, group it by Author, and then use the mean_qi function to calculate these values
forest.data.summary <- group_by(forest.data, author_year) %>% 
  mean_qi(b_Intercept)

# let's plot!
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

# ggsave("results/forestplot_logBLL_no-Halabicky-Iglesias-Min.png", width = 25, height = 15, units = "cm")

# extract draws for EBD assessment, using spread_draws ----
posterior_summary(m.brm)

draws_pooled_b_sd <- spread_draws(m.brm, b_Intercept, sd_author_year__Intercept)


# write.csv(draws_pooled_b_sd, "/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/results/draws_pooled_b_sd_logBLL_no-Halabicky-Iglesias-Min.csv", row.names = F)

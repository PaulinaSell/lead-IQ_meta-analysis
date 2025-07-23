# Bayesian Meta Analysis: Learning
# Harrer, M., Cuijpers, P., Furukawa, T.A., & Ebert, D.D. (2021). Chapter 13: Bayesian Meta-Analysis. In: Doing Meta-Analysis with R: A Hands-On Guide. Boca Raton, FL and London: Chapman & Hall/CRC Press. ISBN 978-0-367-61007-4.

rm(list=ls(all=T))

library(brms)
library(tidyverse)
library(tidybayes)
library(ggridges)
library(glue)
library(stringr)
library(forcats)
library(extraDistr)


load("data/thirdwave.rda")

# In many meta-analyses, τ (the square root of τ2) lies somewhere around 0.3
phcauchy(0.3, sigma = 0.5) # check probability of less than 0.3 between-study heterogeneity τ in half cauchy distribution with sigma 0.5 

priors <- c(prior(normal(0,1), class = Intercept), # overall effect size µ; we grant ~95% prior probability that the true pooled effect size μ lies between −2.0 and 2.0
            prior(cauchy(0,0.5), class = sd)) # between-study heterogeneity τ

m.brm <- brm(TE|se(seTE) ~ 1 + (1|Author),
             data = ThirdWave,
             prior = priors,
             iter = 4000)

# Check Rhat values
summary(m.brm)

# Posterior predictive check
pp_check(m.brm)

# Interpret results
summary(m.brm)
ranef(m.brm)

# lets look at posterior distribution 
post.samples <- posterior_samples(m.brm, c("^b", "^sd"))
names(post.samples)
names(post.samples) <- c("smd", "tau") # rename

# generate density plot for posterior distributions
ggplot(aes(x = smd), data = post.samples) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(post.samples$smd)) +
  labs(x = expression(italic(SMD)),
       y = element_blank()) +
  theme_minimal()

ggplot(aes(x = tau), data = post.samples) +
  geom_density(fill = "lightgreen",               # set the color
               color = "lightgreen", alpha = 0.7) +  
  geom_point(y = 0, 
             x = mean(post.samples$tau)) +        # add point at mean
  labs(x = expression(tau),
       y = element_blank()) +
  theme_minimal()

# check exact probability of effect (SMD) being smaller than certain value (0.3 here)
smd.ecdf <- ecdf(post.samples$smd)
smd.ecdf(0.3)
# 0.002125 -> 0.2%

# lets generate a forest plot!
# first, we have to prepare the data and extract the posterior distribution for each individual study 
study.draws <- spread_draws(m.brm, r_Author[Author,], b_Intercept) %>% 
  mutate(b_Intercept = r_Author + b_Intercept)

# Next: generate distribution of the pooled effect
pooled.effect.draws <- spread_draws(m.brm, b_Intercept) %>% 
  mutate(Author = "Pooled Effect")

# Next: bind study.draws and pooled.effect.draws to one data frame. 
# We then start a pipe again, calling ungroup first, and then use mutate to (1) clean the study labels (i.e. replace dots with spaces), and (2) reorder the study factor levels by effect size (high to low). 
# The result is the data we need for plotting, which we save as forest.data.
forest.data <- bind_rows(study.draws, 
                         pooled.effect.draws) %>% 
  ungroup() %>%
  mutate(Author = str_replace_all(Author, "[.]", " ")) %>% 
  mutate(Author = reorder(Author, b_Intercept))
# Lastly: forest plot should also display the effect size (SMD and credible interval) of each study. 
# We use forest.data data set, group it by Author, and then use the mean_qi function to calculate these values
forest.data.summary <- group_by(forest.data, Author) %>% 
  mean_qi(b_Intercept)

# let's plot!
ggplot(aes(b_Intercept, 
           relevel(Author, "Pooled Effect", 
                   after = Inf)), 
       data = forest.data) +
  
  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = fixef(m.brm)[1, 1], 
             color = "grey", linewidth = 1) +
  geom_vline(xintercept = fixef(m.brm)[1, 3:4], 
             color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "black", 
             linewidth = 1) +
  
  # Add densities
  geom_density_ridges(fill = "lightblue3", 
                      rel_min_height = 0.01, 
                      col = NA, scale = 1,
                      alpha = 0.8) +
  geom_pointinterval(data = forest.data.summary,
                     aes(xmin = .lower, xmax = .upper),
                     linewidth = 1) +
  
  # Add text and labels
  geom_text(data = mutate_if(forest.data.summary, 
                             is.numeric, round, 2),
            aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"), 
                x = Inf), hjust = "inward") +
  labs(x = "Standardized Mean Difference",
       y = element_blank()) +
  theme_minimal()

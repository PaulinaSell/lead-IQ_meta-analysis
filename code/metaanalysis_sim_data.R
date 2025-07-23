# Simulate some data to test model (& compare freq with Bayes?)

rm(list=ls(all=T))

library(brms)
library(tidyverse)
library(tidybayes)
library(ggridges)
library(glue)
library(extraDistr)
library(metafor)

set.seed(177)

author_year <- paste0("author ", sample(2006:2022, 14, replace = F))
beta <- rt(14, 10, -0.5)
se <- rnorm(14, 1, 0.3)
data = data.frame(author_year, beta, se)

ggplot(data, aes(x = factor(author_year, levels = author_year), y = beta)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Set priors
priors <- c(prior(normal(0,1), class = Intercept), # overall effect size µ
            prior(cauchy(0,0.5), class = sd)) # between-study heterogeneity τ

# Prior predictive check
fitPrior <- brm(beta|se(se) ~ 1 + (1|author_year), 
                data = data, 
                prior = priors,
                sample_prior = "only",
                iter = 4000)

# Plot
pp_check(fitPrior)


# Main model
m.brm <- brm(
  beta|se(se) ~ 1 + (1|author_year),
  data = data,
  prior = priors,
  #control =                   # Whenever you see the warning "There were x divergent transitions after warmup." you should really think about increasing adapt_delta.
   # list(adapt_delta = 0.99), # Increasing adapt_delta will slow down the sampler but will decrease the number of divergent transitions threatening the validity of your posterior draws
 # save_pars = save_pars(all = TRUE),
  iter = 4000)

# plot the MCMC chains as well as the posterior distributions
plot(m.brm)

# leave-one-out cross-validation
loo(m.brm, reloo = TRUE) # study 1, 2 and 12 problematic?

# Posterior predictive check
pp_check(m.brm)

# Examine the pairs() plot to diagnose sampling problems
pairs(m.brm)

# Check Rhat values & interpret results
summary(m.brm)
ranef(m.brm)

# lets look at posterior distribution 
str(as_draws_df(m.brm))
post.samples <- as_draws_df(m.brm, variable = c("b_Intercept", "sd_author_year__Intercept"))


# generate density plot for posterior distributions
ggplot(aes(x = b_Intercept), data = post.samples) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +  
  geom_point(y = 0,                               # add point at mean
             x = mean(post.samples$b_Intercept)) +
  labs(x = expression(italic(b_Intercept)),
       y = element_blank()) +
  theme_minimal()

ggplot(aes(x = sd_author_year__Intercept), data = post.samples) +
  geom_density(fill = "lightgreen",               # set the color
               color = "lightgreen", alpha = 0.7) +  
  geom_point(y = 0, 
             x = mean(post.samples$sd_author_year__Intercept)) +  # add point at mean
  labs(x = expression(sd_author_year__Intercept),
       y = element_blank()) +
  theme_minimal()

# now we can check exact probability of effect being smaller than certain value (empirical cumulative distribution function)
smd.ecdf <- ecdf(post.samples$b_Intercept)
smd.ecdf(-0.5)


# lets generate a forest plot!
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
             color = "gray50", linewidth = 0.8) + 
  geom_vline(xintercept = fixef(m.brm)[1, 3:4], # lines for Q2.5 & Q97.5 from pooled effect (95 % CI)
             color = "gray50", linewidth = 0.8, linetype = 2) + 
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
  theme_light()
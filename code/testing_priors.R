# Testing Priors for Bayesian Meta Analysis of Epidemiological Studies on Lead and IQ loss
# run after main-model_log.R

# how to set priors?
# main effect, beta
# visualize
mean <- -1
sd <- 2

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

# heterogeneity, tau with cauchy
sigma = 0.2
phcauchy(0.25, sigma = sigma) # check probability of less than 0.2 between-study heterogeneity τ in half-cauchy distribution with sigma 0.25 
inverseCDF(c(0.025, 0.975), phcauchy, sigma = sigma) # use inverse Cumulative Density Function to find Q2.5 and Q97.5 of the half-cauchy with sigma x

x <- seq(0, 5, by = 0.01)
y <- dhcauchy(x, sigma = sigma)

ggplot() +
  aes(x, y) +
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


# Trying different priors ----
m.brm1 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(normal(-2, 1), class = Intercept), # overall effect size mu
            prior(normal(0, 1), class = sd, lb = 0)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000)

m.brm2 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(normal(-1, 1), class = Intercept), # overall effect size mu
            prior(normal(0, 1), class = sd, lb = 0)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000)

m.brm3 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(normal(-1, 2), class = Intercept), # overall effect size mu
            prior(normal(0, 1), class = sd, lb = 0)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000)

m.brm4 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(normal(-1, 2), class = Intercept), # overall effect size mu
            prior(normal(0, 2), class = sd, lb = 0)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000)

m.brm5 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(uniform(-5, 2), class = Intercept), # overall effect size mu
            prior(uniform(0, 3), class = sd, lb = 0)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000)

m.brm6 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year),
  data = data,
  prior = c(prior(student_t(3, -2.9, 2.8), class = Intercept), # overall effect size mu
            prior(student_t(3, 0, 2.8), class = sd)), #heterogeneity tau
  save_pars = save_pars(all = TRUE),
  iter = 4000)


fitPrior4 <- brm(
  beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), 
  data = data, 
  prior = c(prior(normal(-1, 2), class = Intercept), # overall effect size mu
            prior(normal(0, 2), class = sd, lb = 0)), #heterogeneity tau
  sample_prior = "only",
  iter = 4000)


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


# what is preset prior?
get_prior(beta_ln|se(se_beta_ln) ~ 1 + (1|author_year), data = data)
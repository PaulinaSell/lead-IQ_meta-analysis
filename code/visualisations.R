rm(list=ls(all=T))

library(tidyverse)
library(tidybayes)
library(brms)
library(bayesplot)
library(ggridges)
library(glue)

# Visualizations for Bayesian Meta Analysis of Epidemiological Studies on Lead and IQ loss

# load data & models ----
data <- read.csv("data/study_data_leadIQloss.csv")

# Subsetting full study base for sensitivity analysis: excluding studies with RoB
data_low_medium <- data[!data$author_year %in% c("Crump 2013", "Earl 2016"), ]
data_low <- data[!data$author_year %in% c("Crump 2013", "Earl 2016", "Lucchini 2012", "Lucchini 2019"), ]


read_model <- function(model_name) {
  model = readRDS(paste0("models/", model_name, ".rds"))
  return(model)
}

m.brm_full <- read_model("m.brm_full")
m.brm_low_medium <- read_model("m.brm_low_medium")
m.brm_low <- read_model("m.brm_low")

fitPrior_full <- read_model("fitPrior_full")
fitPrior_low_medium <- read_model("fitPrior_low_medium")
fitPrior_low <- read_model("fitPrior_low")

freq_full <- read_model("freq_full")
freq_low_medium <- read_model("freq_low_medium")
freq_low <- read_model("freq_low")

# create lists ----
# simple list for Bayes MA
result_models <- list(
  full = m.brm_full,
  low_medium = m.brm_low_medium,
  low = m.brm_low
)

# list for "prior, data & posterior" plot
result_configs <- list(
  full = list(
    model = m.brm_full,
    prior = fitPrior_full,
    data = data
  ),
  low_medium = list(
    model = m.brm_low_medium,
    prior = fitPrior_low_medium,
    data = data_low_medium
  ),
  low = list(
    model = m.brm_low,
    prior = fitPrior_low,
    data = data_low
  )
)

# for freq MA
freq_result_models <- list(
  full = freq_full,
  low_medium = freq_low_medium,
  low = freq_low
)

# Priors ----
# beta
x <- seq(-6, 4, by = 0.1)
y <- dnorm(x, -1, 2)
ggplot() +
  aes(x, y) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_area(fill = "steelblue3", alpha = 0.6) +
  theme_minimal() +
  labs(x = NULL, y = NULL)

ggsave(filename = "manuscript/figures/prior_beta.png", 
       width = 20, 
       height = 10, 
       units = "cm")

# tau
a <- seq(0, 10, by = 0.1)
b <- dnorm(a, 0, 2)

ggplot() +
  aes(a, b) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_area(fill = "orange", alpha = 0.6) +
  theme_minimal() +
  labs(x = NULL, y = NULL)

ggsave(filename = "manuscript/figures/prior_tau.png", 
       width = 20, 
       height = 10, 
       units = "cm")

# Traceplot incl. warmup (no loop) ----
 posterior_samples_warm = as_draws_df(m.brm_full, inc_warmup = T)
 names(posterior_samples_warm)[names(posterior_samples_warm) == "b_Intercept"] = "beta"
 names(posterior_samples_warm)[names(posterior_samples_warm) == "sd_author_year__Intercept"] = "tau"

 mcmc_trace(posterior_samples_warm,
            pars = c("beta", "tau"),
            facet_args = list(ncol = 1)) +
   vline_at(2000, color = "red", linetype = 2, size = 0.5)  # mark end of warmup
 
 

# Posterior: beta & tau ----

# loop -> plot & save posterior for all 3 data configs for beta & tau
for (model_name in names(result_models)) {
   model <- result_models[[model_name]]
 
 post_samples <- as_draws_df(model, variable = c("b_Intercept", "sd_author_year__Intercept"))
 
 beta <- ggplot(aes(x = b_Intercept), data = post_samples) +
   geom_density(fill = "steelblue",
                color = "steelblue", alpha = 0.7) +  
   # geom_vline(xintercept = mean(post_samples$b_Intercept), 
   # linetype = "dotted", 
   # color = "red") +
   # geom_vline(xintercept = Mode(post_samples$b_Intercept), 
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
 
 tau <- ggplot(aes(x = sd_author_year__Intercept), data = post_samples) +
   geom_density(fill = "lightgreen",
                color = "lightgreen", alpha = 0.7) +  
   # geom_vline(xintercept = mean(post_samples$sd_author_year__Intercept), 
   #            linetype = "dotted", 
   #            color = "red") +
   # geom_vline(xintercept = Mode(post_samples$sd_author_year__Intercept), 
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
 
 ggsave(filename = paste0("results/posterior_b_", model_name, ".png"), 
        plot = beta, 
        width = 30, 
        height = 10, 
        units = "cm")
 
 ggsave(filename = paste0("results/posterior_tau_", model_name, ".png"), 
        plot = tau, 
        width = 30, 
        height = 10, 
        units = "cm")

}
 
 # Prior, data & posterior in one plot: looped over all 3 sensitivity analyses----

 for (config_name in names(result_configs)) {
   current_config <- result_configs[[config_name]]
   model <- current_config$model
   prior <- current_config$prior
   data <- current_config$data
   

 # Beta
 # Prepare data
 posterior_samples <- as_draws_df(model)
 prior_samples <- as_draws_df(prior)
 
 plot_data_combined <- data.frame(
   value = c(prior_samples[["b_Intercept"]], 
             posterior_samples[["b_Intercept"]],
             data$beta_ln),
   distribution = factor(c(rep("Prior", nrow(prior_samples)),
                           rep("Posterior", nrow(posterior_samples)),
                           rep("Observed Effects", nrow(data))
                           ),
                         levels = c("Prior", "Observed Effects", "Posterior"))
 )
 
 beta <- 
   ggplot(plot_data_combined, aes(x = value, fill = distribution)) +
   geom_density(alpha = 0.5, linewidth = 0.3) +
   labs(
     # title = "Prior, Observed Effects, and Posterior of Beta",
     x = "Beta", 
     y = "Density",
     fill = "Distribution") +
   scale_fill_manual(values = c("Prior" = "#e9c46a", "Posterior" = "#2a9d8f", "Observed Effects" = "#e76f51"),
                     name = "") +
   theme(
     legend.position.inside = c(.95, .95),
     legend.justification = c("right", "top"),
     panel.border = element_blank(),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     panel.background = element_blank()
   ) 
 
 # Heterogeneity:

 plot_data_combined <- data.frame(
   value = c(prior_samples[["sd_author_year__Intercept"]], 
             posterior_samples[["sd_author_year__Intercept"]]),
   distribution = factor(c(rep("Prior", nrow(prior_samples)),
                           rep("Posterior", nrow(posterior_samples))),
                         levels = c("Prior", "Posterior"))
 )
 
 heterogeneity <- 
   ggplot(plot_data_combined, aes(x = value, fill = distribution)) +
   geom_density(alpha = 0.5, linewidth = 0.3) +
   labs(
     # title = "Prior, Observed Effects, and Posterior of Beta",
     x = "Tau", 
     y = "Density",
     fill = "Distribution") +
   scale_fill_manual(values = c("Prior" = "#e9c46a", "Posterior" = "#2a9d8f"),
                     name = "") +
   theme(
     legend.position = c(.95, .95),
     legend.justification = c("right", "top"),
     panel.border = element_blank(),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     panel.background = element_blank()
   ) 
 
 plot_list <- list(
   heterogeneity = heterogeneity,
   beta = beta
 )
 
 for(p_name in names(plot_list)) {
   
 ggsave(
   filename = paste0("results/prior_obs-effects_posterior_", p_name, "-", config_name, ".png"),
   plot = plot_list[[p_name]],
   width = 22,
   height = 15,
   units = "cm")
 }
 }

 # Freq MA forest plot looping ----
 for (freq_model_name in names(freq_result_models)) {
   model <- freq_result_models[[freq_model_name]]
   
   # plot & save forests
   png(paste0("results/freq_forestplot_", freq_model_name, ".png"), 
       width = 20, height = 15, units = "cm", res = 300)
   
   forest(model, 
          slab = model$data$author_year,
          main = "Random Effects Meta-Analysis")
   
   dev.off()
  }

 
 # Bayesian MA forest plot looping ----
 
 for (model_name in names(result_models)) {
   model <- result_models[[model_name]] 
 
 # 1. extract posterior distribution for each individual study 
 study.draws <- spread_draws(model, r_author_year[author_year,], b_Intercept) %>% 
   mutate(b_Intercept = r_author_year + b_Intercept)
 
 # 2. generate distribution of the pooled effect
 pooled.effect.draws <- spread_draws(model, b_Intercept) %>% 
   mutate(author_year = "Pooled Effect")
 
 # 3. bind study.draws and pooled.effect.draws to one data frame. 
 forest.data <- bind_rows(study.draws, 
                          pooled.effect.draws) %>% 
   ungroup() %>%
   mutate(author_year = str_replace_all(author_year, "[.]", " ")) %>%
   mutate(author_year = fct_rev(factor(author_year)))
 
 # 4. forest plot should also display the effect size (SMD and credible interval) of each study. 
 # We use forest.data data set, group it by author, and then use the mean_qi function to calculate these values
 forest.data.summary <- group_by(forest.data, author_year) %>% 
   mean_qi(b_Intercept)
 
 p <- ggplot(aes(b_Intercept, 
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
   geom_vline(xintercept = fixef(model)[1, 1], # line at the pooled effect estimate
              color = "gray50", linewidth = 0.8, linetype = 2) + 
   geom_vline(xintercept = 0, color = "gray20", 
              linewidth = 1) + # line at zero (null effect line)
   
   # Add density ridges
   geom_density_ridges(fill = "steelblue3", # ridge density plots for each study/row
                       rel_min_height = 0.01, # removes very small density values
                       col = NA, scale = 1, # no outline color, scale of distributions
                       alpha = 0.8) +
   
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
   theme_light() +
   theme(panel.border = element_blank())
 
 ggsave(
   filename = paste0("results/forestplot_", model_name, ".png"),
   plot = p,
   width = 25,
   height = 15, 
   units = "cm")
 
 }
 
 
# Plotting curves: lead & IQ loss function with uncertainty ----
 
# Bayes:
 for (model_name in names(result_models)) {
   model <- result_models[[model_name]] 
   
 # extract results from brmsfit object
 mean_beta <- posterior_summary(model, variable = "b_Intercept")[, "Estimate"] 
 lower_beta <- posterior_summary(model, variable = "b_Intercept")[, 3]
 upper_beta <- posterior_summary(model, variable = "b_Intercept")[, 4]
 
 IQ_loss <- data.frame(
   blood_lead = seq(0, 30, by = 0.1)
 )
 
 IQ_loss$mean <- mean_beta * log(IQ_loss$blood_lead + 1)
 IQ_loss$lowerCI <- lower_beta * log(IQ_loss$blood_lead + 1)
 IQ_loss$upperCI <- upper_beta * log(IQ_loss$blood_lead + 1)
 
 plot <- ggplot(IQ_loss, aes(x = blood_lead)) +
   geom_line(aes(y = mean, color = "Mean"), show.legend = F) +
   geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), fill = "lightblue", alpha = 0.4) +
   labs(
     x = "Blood Lead Level (µg/dL)",
     y = "IQ Loss (points)"
   ) +
   xlim(0, 30) +
   ylim(-14, 0) +
   theme_minimal() 
 
 ggsave(
   filename = paste0("results/Bayes_curve_", model_name, ".png"),
   plot = plot,
   width = 22,
   height = 15, 
   units = "cm")
 }
 
 # Freq: 
 for (freq_model_name in names(freq_result_models)) {
   model <- freq_result_models[[freq_model_name]]
   
   # extract results from brmsfit object
   mean_beta <- model$beta[1, 1]
   lower_beta <- model$ci.lb[1]
   upper_beta <- model$ci.ub[1]
   
   IQ_loss <- data.frame(
     blood_lead = seq(0, 30, by = 0.1)
   )
   
   IQ_loss$mean <- mean_beta * log(IQ_loss$blood_lead + 1)
   IQ_loss$lowerCI <- lower_beta * log(IQ_loss$blood_lead + 1)
   IQ_loss$upperCI <- upper_beta * log(IQ_loss$blood_lead + 1)
   
   plot <- ggplot(IQ_loss, aes(x = blood_lead)) +
     geom_line(aes(y = mean, color = "Mean"), show.legend = F) +
     geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), fill = "lightblue", alpha = 0.4) +
     labs(
       x = "Blood Lead Level (µg/dL)",
       y = "IQ Loss (points)"
     ) +
     xlim(0, 30) +
     ylim(-14, 0) +
     theme_minimal() 
   
   ggsave(
     filename = paste0("results/freq_curve_", freq_model_name, ".png"),
     plot = plot,
     width = 22,
     height = 15, 
     units = "cm")
   
 }
 
 
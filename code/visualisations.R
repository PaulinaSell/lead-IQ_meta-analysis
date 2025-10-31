# Visualizations for Bayesian Meta Analysis of Epidemiological Studies on Lead and IQ loss
# run after main-model_log.R

# looking at the data 
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
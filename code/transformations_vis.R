# Compare log vs linear Exposure-Response Functions ----
# (both are beta from Crump et al. 2013, original is fpr ln(BLL), beta_lin was transformed to linear)

library(ggplot2)

beta_log <- -2.650
beta_linear <- -0.277

blood_lead <- seq(0, 30, by = 0.1)

IQ_loss_log <- beta_log * log(blood_lead + 1)
IQ_loss_linear <- beta_linear * blood_lead

plot_data <- data.frame(
  blood_lead = blood_lead,
  log_ERF = IQ_loss_log,
  linear_ERF = IQ_loss_linear
)

ggplot(plot_data, aes(x = blood_lead)) +
  geom_line(aes(y = log_ERF, color = "Logarithmic")) +
  geom_line(aes(y = linear_ERF, color = "Linear")) +
  scale_color_manual(values = c("Logarithmic" = "blue", "Linear" = "red")) +
  labs(
    x = "Blood Lead Level (µg/dL)",
    y = "IQ Loss (points)"
  ) +
  theme_minimal()

# transformation to linear appears to be invalid

# visually check transformations to ln(BLL) as an alternative
# case 1: from log10(BLL) to ln(BLL): Dantzer et al. 2020
Dantzerbeta_log10 <- -13.15
Dantzerbeta_ln <- -5.711

blood_lead <- seq(0, 30, by = 0.1)

IQ_loss_log10 <- Dantzerbeta_log10 * log10(blood_lead + 1)
IQ_loss_ln <- Dantzerbeta_ln * log(blood_lead + 1)

Dantzerplot_data <- data.frame(
  blood_lead = blood_lead,
  log_ERF = IQ_loss_log10,
  ln_ERF = IQ_loss_ln
)

ggplot(Dantzerplot_data, aes(x = blood_lead)) +
  geom_line(aes(y = log_ERF, color = "log(10)")) +
  geom_line(aes(y = ln_ERF, color = "ln")) +
  scale_color_manual(values = c("log(10)" = "blue", "ln" = "red")) +
  labs(
    x = "Blood Lead Level (µg/dL)",
    y = "IQ Loss (points)"
  ) +
  theme_minimal()
# very similar

# case 2: from log2(BLL) to ln(BLL): Desrochers-Couture et al. 2018
Desrochersbeta_log2 <- 0.014
Desrochersbeta_ln <- 0.020

blood_lead <- seq(0, 30, by = 0.1)

IQ_loss_log2 <- Desrochersbeta_log2 * log2(blood_lead + 1)
IQ_loss_ln <- Desrochersbeta_ln * log(blood_lead + 1)

Desrochersplot_data <- data.frame(
  blood_lead = blood_lead,
  log_ERF = IQ_loss_log2,
  ln_ERF = IQ_loss_ln
)

ggplot(Desrochersplot_data, aes(x = blood_lead)) +
  geom_line(aes(y = log_ERF, color = "log(2)")) +
  geom_line(aes(y = ln_ERF, color = "ln")) +
  scale_color_manual(values = c("log(2)" = "blue", "ln" = "red")) +
  labs(
    x = "Blood Lead Level (µg/dL)",
    y = "IQ Loss (points)"
  ) +
  theme_minimal()
# very similar

# case 3: from linear to ln(BLL): Iglesias et al. 2011, reference BLL 2.2 µg/dL
Iglesiasbeta_linear <- -0.940
Iglesiasbeta_ln <- -2.068

blood_lead <- seq(0, 30, by = 0.1)

IQ_loss_linear <- Iglesiasbeta_linear * blood_lead
IQ_loss_ln <- Iglesiasbeta_ln * log(blood_lead + 1)

Iglesiasplot_data <- data.frame(
  blood_lead = blood_lead,
  linear_ERF = IQ_loss_linear,
  ln_ERF = IQ_loss_ln
)

ggplot(Iglesiasplot_data, aes(x = blood_lead)) +
  geom_line(aes(y = linear_ERF, color = "linear")) +
  geom_line(aes(y = ln_ERF, color = "ln")) +
  scale_color_manual(values = c("linear" = "blue", "ln" = "red")) +
  labs(
    x = "Blood Lead Level (µg/dL)",
    y = "IQ Loss (points)"
  ) +
  theme_minimal()
# does not appear to be a valid approximation in high concentration, lets zoom in on lower concentration

blood_lead_low <- seq(0, 5, by = 0.1)

IQ_loss_linear <- Iglesiasbeta_linear * blood_lead_low
IQ_loss_ln <- Iglesiasbeta_ln * log(blood_lead_low + 1)

Iglesiasplot_data <- data.frame(
  blood_lead_low = blood_lead_low,
  linear_ERF = IQ_loss_linear,
  ln_ERF = IQ_loss_ln
)

ggplot(Iglesiasplot_data, aes(x = blood_lead_low)) +
  geom_line(aes(y = linear_ERF, color = "linear")) +
  geom_line(aes(y = ln_ERF, color = "ln")) +
  scale_color_manual(values = c("linear" = "blue", "ln" = "red")) +
  labs(
    x = "Blood Lead Level (µg/dL)",
    y = "IQ Loss (points)"
  ) +
  theme_minimal()
# does not appear to be a valid approximation in lower concentration

# Faceted prior/posterior figures for the prior sensitivity analyses.
# Run from the repository root: Rscript code/figures_priors_oneway.R

library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(posterior)

source("code/prior_configs.R")

FIGURE_DIR <- "results/sensitivity/priors"
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

data_full <- read.csv("data/study_data_leadIQloss.csv")

DIST_COLOURS <- c(
  "Prior" = "#e9c46a",
  "Observed effects" = "#e76f51",
  "Posterior" = "#2a9d8f"
)

config_levels <- function(configs) map_chr(configs, "label")

mu_draws <- function(configs) {
  map_dfr(configs, function(cfg) {
    prior_samples <- as_draws_df(load_prior(cfg))
    posterior_samples <- as_draws_df(load_posterior(cfg))

    bind_rows(
      tibble(value = prior_samples$b_Intercept, distribution = "Prior"),
      tibble(value = data_full$beta_ln, distribution = "Observed effects"),
      tibble(value = posterior_samples$b_Intercept, distribution = "Posterior")
    ) %>%
      mutate(config = cfg$label)
  }) %>%
    mutate(
      config = factor(config, levels = config_levels(configs)),
      distribution = factor(distribution, levels = names(DIST_COLOURS))
    )
}

tau_draws <- function(configs) {
  map_dfr(configs, function(cfg) {
    prior_samples <- as_draws_df(load_prior(cfg))
    posterior_samples <- as_draws_df(load_posterior(cfg))

    bind_rows(
      tibble(
        value = prior_samples$sd_author_year__Intercept,
        distribution = "Prior"
      ),
      tibble(
        value = posterior_samples$sd_author_year__Intercept,
        distribution = "Posterior"
      )
    ) %>%
      mutate(config = cfg$label)
  }) %>%
    mutate(
      config = factor(config, levels = config_levels(configs)),
      distribution = factor(distribution, levels = names(DIST_COLOURS))
    )
}

base_theme <- function() {
  theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(hjust = 0, face = "bold")
    )
}

plot_mu <- function(draws, xlim = c(-8, 5)) {
  ggplot(draws, aes(x = value, fill = distribution)) +
    geom_density(alpha = 0.5, linewidth = 0.3) +
    facet_wrap(~config, ncol = 1) +
    coord_cartesian(xlim = xlim) +
    scale_fill_manual(values = DIST_COLOURS, name = NULL) +
    labs(x = "Pooled effect, mu", y = "Density") +
    base_theme()
}

plot_mu_and_tau <- function(mu, tau) {
  combined <- bind_rows(
    mutate(mu, parameter = "Pooled effect, mu"),
    mutate(tau, parameter = "Heterogeneity, tau")
  ) %>%
    mutate(
      parameter = factor(
        parameter,
        levels = c("Pooled effect, mu", "Heterogeneity, tau")
      )
    )

  ggplot(combined, aes(x = value, fill = distribution)) +
    geom_density(alpha = 0.5, linewidth = 0.3) +
    facet_wrap(~ config + parameter, ncol = 2, scales = "free") +
    scale_fill_manual(values = DIST_COLOURS, name = NULL) +
    labs(x = NULL, y = "Density") +
    base_theme()
}

save_svg <- function(plot, filename, height) {
  ggsave(
    filename = file.path(FIGURE_DIR, filename),
    plot = plot,
    width = 18,
    height = height,
    units = "cm"
  )
}

effect_configs <- keep(oneway_prior_configs, ~ .x$arm == "Effect prior")
heterogeneity_configs <- keep(
  oneway_prior_configs,
  ~ .x$arm == "Heterogeneity prior"
)

save_svg(plot_mu(mu_draws(joint_prior_configs)), "mu_joint-priors.svg", 18)
save_svg(plot_mu(mu_draws(effect_configs)), "mu_effect-prior.svg", 18)
save_svg(
  plot_mu_and_tau(
    mu_draws(heterogeneity_configs),
    tau_draws(heterogeneity_configs)
  ),
  "mu-tau_heterogeneity-prior.svg",
  20
)

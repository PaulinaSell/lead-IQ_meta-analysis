# forest plot -> make similar to metafor

study.draws <- spread_draws(
  m.brm,
  r_author_year[author_year, ],
  b_Intercept
) %>%
  mutate(b_Intercept = r_author_year + b_Intercept)

# Generate distribution of the pooled effect
pooled.effect.draws <- spread_draws(m.brm, b_Intercept) %>%
  mutate(author_year = "Random-Effects Model")

forest.data <- bind_rows(study.draws, pooled.effect.draws) %>%
  ungroup() %>%
  mutate(author_year = str_replace_all(author_year, "[.]", " ")) %>%
  mutate(author_year = fct_rev(factor(author_year)))

forest.data.summary <- group_by(forest.data, author_year) %>%
  mean_qi(b_Intercept)


# let's plot!
ggplot(
  aes(b_Intercept, relevel(author_year, "Random-Effects Model", after = Inf)),
  data = forest.data
) +
  geom_vline(xintercept = 0, linewidth = 0.6, linetype = 2) +
  geom_density_ridges(
    fill = "steelblue3",
    rel_min_height = 0.01,
    col = NA,
    scale = 1,
    alpha = 0.8
  ) +
  geom_pointinterval(
    data = forest.data.summary,
    aes(xmin = .lower, xmax = .upper),
    fatten_point = 1.5,
    linewidth = 0.8
  ) +
  geom_text(
    data = mutate_if(
      forest.data.summary,
      is.numeric,
      round,
      2
    ),
    aes(label = glue("{format(b_Intercept, nsmall = 2)} [{format(.lower, nsmall = 2)}, {format(.upper, nsmall = 2)}]"), x = 3.5)
    hjust = 0.5
  ) +
  labs(
    title = "Random Effects Meta-Analysis",
    x = "Observed Outcome",
    y = element_blank()
  ) +
  xlim(-15, 5) +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(hjust = 0)
  )

# forest plot -> make similar to metafor

study.draws <- spread_draws(
  m.brm,
  r_author_year[author_year, ],
  b_Intercept
) %>%
  mutate(b_Intercept = r_author_year + b_Intercept)

# Generate distribution of the pooled effect
pooled.effect.draws <- spread_draws(m.brm, b_Intercept) %>%
  mutate(author_year = "Pooled Estimate")

forest.data <- bind_rows(study.draws, pooled.effect.draws) %>%
  ungroup() %>%
  mutate(author_year = str_replace_all(author_year, "[.]", " ")) %>%
  mutate(author_year = fct_rev(factor(author_year)))

forest.data.summary <- group_by(forest.data, author_year) %>%
  mean_qi(b_Intercept)

# ### CHANGE 1: compute how many levels the y-axis has, needed for annotation positioning ###
n_studies <- nlevels(forest.data$author_year) # includes "Pooled Estimate"


# let's plot!
ggplot(
  aes(b_Intercept, relevel(author_year, "Pooled Estimate", after = Inf)),
  data = forest.data
) +
  geom_vline(xintercept = 0, linewidth = 0.6, linetype = "dotted") +
  geom_density_ridges(
    fill = "grey",
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
    aes(
      label = glue(
        "{format(b_Intercept, nsmall = 2)} [{format(.lower, nsmall = 2)}, {format(.upper, nsmall = 2)}]"
      ),
      x = 3.5
    ),
    hjust = 0.2
  ) +
  # ### CHANGE 2: subtitle-style column headers drawn as plot annotations ###
  annotate(
    "text",
    x = -15,
    y = n_studies + 1, # left-aligned, just above the top horizontal line
    label = "Study",
    hjust = 0,
    fontface = "bold",
    size = 3.5
  ) +
  annotate(
    "text",
    x = 3.5,
    y = n_studies + 1, # right side, matching the estimates column
    label = "Estimate [95% CrI]",
    hjust = 0.2,
    fontface = "bold",
    size = 3.5
  ) +
  # ### CHANGE 3: horizontal lines mimicking metafor style ###
  # Line below the column headers / above the first study row
  annotate(
    "segment",
    x = -15,
    xend = 5,
    y = n_studies + 0.5,
    yend = n_studies + 0.5,
    linewidth = 0.5
  ) +
  # Line above the pooled estimate row
  annotate(
    "segment",
    x = -15,
    xend = 5,
    y = 1.5,
    yend = 1.5,
    linewidth = 0.5
  ) +
  labs(
    title = "Bayesian Random Effects Meta-Analysis",
    x = "IQ shift per log-unit increase in BLL (µg/dL)",
    y = element_blank()
  ) +
  scale_x_continuous(
    limits = c(-15, 5),
    breaks = seq(-15, 5, by = 5), # explicit tick marks
    expand = expansion(add = c(0, 4)) # 0 padding left, 2 units right
  ) +
  # ### CHANGE 4: y-axis expanded upward to give room for the header annotations ###
  scale_y_discrete(
    expand = expansion(add = c(0.5, 2)) # 0.5 padding bottom, 2 units top
  ) +
  guides(x = guide_axis(cap = "both")) + # cap both ends
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(t = 40, r = 10, b = 40, l = 20), # increase top margin
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(hjust = 0, margin = margin(r = -70)), # negative right margin pulls text right
    axis.line.x = element_line(linewidth = 0.5), # draw the x axis line
    axis.ticks.x = element_line(linewidth = 0.5), # draw tick marks
    axis.ticks.length.x = unit(0.3, "cm"),
  )

forest(
  freq_model,
  slab = freq_model$data$author_year,
  mlab = "Pooled Estimate",
  xlab = "IQ shift per log-unit increase in BLL (µg/dL)",
  main = "Frequentist Random Effects Meta-Analysis"
)

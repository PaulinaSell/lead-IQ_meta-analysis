# forest plot -> make similar to metafor

library(posterior)
library(HDInterval)
library(brms)
library(metafor)
library(tidyverse)
library(tidybayes)
library(ggridges)
library(glue)

m.brm <- readRDS(
  "/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/models/main/m.brm_full.rds"
)
fitPrior <- readRDS(
  "/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/models/main/fitPrior_full.rds"
)
freq_model <- readRDS(
  "/Users/paulinasell/Documents/UBA/PARC/Metaanalysis_lead_IQloss/RProj/models/main/freq_full.rds"
)


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

# Compute how many levels the y-axis has, needed for annotation positioning
n_studies <- nlevels(forest.data$author_year) # includes "Pooled Estimate"


# let's plot!
Bayes_forest_plot <- ggplot(
  aes(b_Intercept, relevel(author_year, "Pooled Estimate", after = Inf)),
  data = forest.data
) +
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = 0,
    yend = n_studies + 0.5,
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  # geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed") +
  geom_density_ridges(
    fill = "grey",
    rel_min_height = 0.01,
    col = NA,
    scale = 1
  ) +
  geom_pointinterval(
    data = forest.data.summary,
    aes(xmin = .lower, xmax = .upper),
    shape = 22,
    point_fill = "black",
    fatten_point = 1.5,
    linewidth = 0.15
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
    hjust = 0.2,
    size = 4
  ) +
  # Subtitle column headers
  annotate(
    "text",
    x = -17.4,
    y = n_studies + 1, # left-aligned, just above the top horizontal line
    label = "Study",
    hjust = 0.6,
    fontface = "bold",
    size = 4.1
  ) +
  annotate(
    "text",
    x = 3.5,
    y = n_studies + 1, # right side, matching the estimates column
    label = "Estimate [95% CrI]",
    hjust = 0.27,
    fontface = "bold",
    size = 4.1
  ) +
  # Horizontal lines mimicking metafor style
  # Line below the column headers / above the first study row & above pooled estimate
  annotate(
    "segment",
    x = -18.5,
    xend = 7.5,
    y = n_studies + 0.5,
    yend = n_studies + 0.5,
    linewidth = 0.25
  ) +
  annotate(
    "segment",
    x = -18.5,
    xend = 7.5,
    y = 1.5,
    yend = 1.5,
    linewidth = 0.25
  ) +
  labs(
    title = "Bayesian Random Effects Meta-Analysis",
    x = "IQ shift per log-unit increase in BLL (µg/dL)",
    y = element_blank()
  ) +
  # y-axis expanded upward to give room for the header annotations
  scale_y_discrete(
    expand = expansion(add = c(0.5, 2.4)) # 0.5 padding bottom, 2 units top
  ) +
  scale_x_continuous(
    # remove limits here — let coord_cartesian control clipping instead
    breaks = seq(-15, 5, by = 5),
    expand = expansion(add = c(0, 4))
  ) +
  coord_cartesian(
    xlim = c(-15, 5),
    clip = "off" # allows segments/annotations to extend into the margin
  ) +
  guides(x = guide_axis(cap = "both")) + # cap both ends
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(t = 20, r = 10, b = 40, l = 20), # increase top margin
    plot.title = element_text(hjust = 0.3, face = "bold", size = 13), # increased title font size
    axis.text.y = element_text(
      hjust = 0,
      margin = margin(r = -70),
      size = 11.5
    ), # increased y-axis label size
    axis.text.x = element_text(vjust = -1, size = 12), # increased x-axis tick label size
    axis.title.x = element_text(vjust = -2.7, size = 11.5),
    axis.line.x = element_line(linewidth = 0.25), # draw the x axis line
    axis.ticks.x = element_line(linewidth = 0.25), # draw tick marks
    axis.ticks.length.x = unit(0.3, "cm"),
  )

ggsave(
  "manuscript/tables_figures/main/Bayes_forest_plot.png",
  Bayes_forest_plot,
  width = 20,
  height = 15,
  units = "cm"
)

png(
  "manuscript/tables_figures/main/freq_forestplot.png",
  width = 20,
  height = 15,
  units = "cm",
  res = 300
)

forest(
  freq_model,
  slab = freq_model$data$author_year,
  mlab = "Pooled Estimate",
  xlab = "IQ shift per log-unit increase in BLL (µg/dL)",
  main = "Frequentist Random Effects Meta-Analysis"
)

dev.off()

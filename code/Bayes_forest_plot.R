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

# Build an explicit factor level order that inserts two blank "spacer" levels
# This replicates the one-line visual gap seen in metafor forest plots
# We derive the study levels from the fct_rev order so their sequence is preserved
study_levels <- setdiff(levels(forest.data$author_year), "Pooled Estimate")

new_levels <- c(
  "Pooled Estimate", # y = 1  (bottom)
  "spacer_bottom", # y = 2  (blank gap)
  study_levels, # y = 3 … n_studies + 1
  "spacer_top" # y = n_studies + 2  (blank gap, just below header)
)

# Apply the new level order to both data frames
forest.data <- forest.data %>%
  mutate(author_year = factor(author_year, levels = new_levels))

forest.data.summary <- forest.data.summary %>%
  mutate(author_year = factor(as.character(author_year), levels = new_levels))

# Append two spacer rows (all-NA estimates) to the summary data.
# geom_pointinterval and geom_text will be filtered to skip these rows, so nothing is drawn at the spacer positions
spacer_rows <- tibble(
  author_year = factor(c("spacer_bottom", "spacer_top"), levels = new_levels),
  b_Intercept = NA_real_,
  .lower = NA_real_,
  .upper = NA_real_
)
forest.data.summary <- bind_rows(forest.data.summary, spacer_rows)

# Total number of discrete y levels (original studies + 2 spacers).
# Use this instead of n_studies for all annotation y-positions below.
n_levels_total <- n_studies + 2

# let's plot!
Bayes_forest_plot <- ggplot(
  aes(b_Intercept, author_year),
  data = forest.data
) +
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = 0,
    yend = n_levels_total,
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  geom_density_ridges(
    fill = "grey",
    rel_min_height = 0.01,
    col = NA,
    scale = 1
  ) +
  geom_pointinterval(
    data = forest.data.summary %>% filter(!is.na(b_Intercept)),
    aes(xmin = .lower, xmax = .upper),
    shape = 22,
    point_fill = "black",
    fatten_point = 1.5,
    linewidth = 0.15
  ) +
  geom_text(
    data = mutate_if(
      forest.data.summary %>% filter(!is.na(b_Intercept)),
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
  # Column headers ("subtitle")
  annotate(
    "text",
    x = -17.4,
    y = n_levels_total + 1,
    label = "Study",
    hjust = 0.6,
    fontface = "bold",
    size = 4.1
  ) +
  annotate(
    "text",
    x = 3.5,
    y = n_levels_total + 1,
    label = "Estimate [95% CrI]",
    hjust = 0.27,
    fontface = "bold",
    size = 4.1
  ) +
  annotate(
    "segment",
    x = -18.5,
    xend = 7.5,
    y = n_levels_total,
    yend = n_levels_total,
    linewidth = 0.25
  ) +
  annotate(
    "segment",
    x = -18.5,
    xend = 7.5,
    y = 2,
    yend = 2,
    linewidth = 0.25
  ) +
  labs(
    title = "Bayesian Random Effects Meta-Analysis",
    x = "IQ shift per log-unit increase in BLL (µg/dL)",
    y = element_blank()
  ) +
  # y-axis expanded upward to give room for the header annotations
  scale_y_discrete(
    expand = expansion(add = c(0.5, 2.4)), # 0.5 padding bottom, 2 units top
    labels = function(x) ifelse(x %in% c("spacer_bottom", "spacer_top"), "", x),
    drop = FALSE
  ) +
  scale_x_continuous(
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
    ),
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

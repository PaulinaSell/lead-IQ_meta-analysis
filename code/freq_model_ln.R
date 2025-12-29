# Frequentist Random Effects Meta Analysis Lead - IQ loss in Children

rm(list=ls(all=T))

library(metafor)

data_full = read.csv("data/study_data_leadIQloss.csv")

# Subsetting full study base for sensitivity analysis: excluding studies with RoB
data_low_medium <- data_full[!data_full$author_year %in% c("Crump 2013", "Earl 2016"), ]
data_low <- data_full[!data_full$author_year %in% c("Crump 2013", "Earl 2016", "Lucchini 2012", "Lucchini 2019"), ]

models_data <- list(
  full = data_full,
  low_medium = data_low_medium,
  low = data_low
)
  
for (data_subset in names(models_data)) {
  data <- models_data[[data_subset]]
  
  # Run model equivalent to Bayesian model (same data, random effects)
  rma_model <- rma(yi = beta_ln, 
                   sei = se_beta_ln,
                   data = data,
                   method = "REML")  # Restricted Maximum Likelihood

  cat("\n===========================================\n")
  cat("Model including:", data_subset, "RoB studies \n")
  cat("===========================================\n")
  print(summary(rma_model))

  png(paste0("results/freq_forestplot_lnBLL_", data_subset, ".png"), 
      width = 20, height = 15, units = "cm", res = 300)

  forest(rma_model, 
        slab = data$author_year,
        main = "Random Effects Meta-Analysis")

  dev.off() # close device

}


# Funnel plot
funnel(rma_model)

# Diagnostics
# Influence diagnostics (equivalent to Pareto-k checks)
inf <- influence(rma_model)
plot(inf)

# Residuals
plot(rma_model)

# Leave-One-Out Analysis
loo_results <- leave1out(rma_model)
print(loo_results)


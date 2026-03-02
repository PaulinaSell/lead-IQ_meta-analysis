# Frequentist Random Effects Meta Analysis Lead - IQ loss in Children

library(metafor)

data_full = read.csv("data/study_data_leadIQloss.csv")

# run model
rma_model <- rma(
  yi = beta_ln,
  sei = se_beta_ln,
  data = data_full,
  method = "REML"
) # Restricted Maximum Likelihood

saveRDS(rma_model, file = "models/main/freq_full.rds")

rma_model <- readRDS("models/main/freq_full.rds")

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

forest(rma_model)

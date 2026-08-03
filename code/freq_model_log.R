# Frequentist Random Effects Meta Analysis Lead - IQ loss in Children

library(metafor)

data_full = read.csv("data/study_data_leadIQloss.csv")

# run unadjusted model
rma_model <- rma(
  yi = beta_ln,
  sei = se_beta_ln,
  data = data_full,
  method = "REML"
) # Restricted Maximum Likelihood

print(rma_model)

rma_model_knha <- rma(
  yi = beta_ln,
  sei = se_beta_ln,
  data = data_full,
  method = "REML",
  test = "knha"
) # Small sample adjustment

print(rma_model_knha)

saveRDS(rma_model_knha, file = "models/main/freq_full.rds")

# Funnel plot
funnel(rma_model_knha)

# Diagnostics
# Influence diagnostics (equivalent to Pareto-k checks)
inf <- influence(rma_model_knha)
plot(inf)

# Residuals
plot(rma_model_knha)

# Leave-One-Out Analysis
# loo_results <- leave1out(rma_model_knha)
# print(loo_results)

forest(rma_model_knha)

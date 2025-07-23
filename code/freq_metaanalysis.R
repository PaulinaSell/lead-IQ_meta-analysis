# Frequentist Random Effects Meta Analysis Lead - IQ loss in Children

rm(list=ls(all=T))

library(metafor)

data = read.csv("data/study_data_leadIQloss.csv")


# Equivalent to bayesian model with random effects

rma_model <- rma(yi = beta_lin, 
                 sei = se_beta_lin, 
                 data = data,
                 method = "REML")  # Restricted Maximum Likelihood

summary(rma_model)


png("results/freq_forestplot1.png", 
    width = 20, height = 15, units = "cm", res = 300)

# Jetzt den Plot erstellen
forest(rma_model, 
       slab = data$author_year,
       main = "Random Effects Meta-Analysis")

# Datei schließen und speichern
dev.off()

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

# Forest plot mit Leave-One-Out results
forest(loo_results)

# Test für Heterogenität
# Q-Test und I²
print(paste("Q =", round(rma_model$QE, 2), 
            "df =", rma_model$k - rma_model$p, 
            "p =", round(rma_model$QEp, 4)))
print(paste("I² =", round(rma_model$I2, 1), "%"))
print(paste("τ² =", round(rma_model$tau2, 4)))

# Alternative: Fixed Effects Model zum Vergleich
rma_fixed <- rma(yi = beta_lin, 
                 sei = se_beta_lin, 
                 data = data,
                 method = "FE")
summary(rma_fixed)

# Modellvergleich
anova(rma_fixed, rma_model)

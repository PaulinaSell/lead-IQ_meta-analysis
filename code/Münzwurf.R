# Bayessche Posterior-Berechnung für Münzwurf
# Prior: N(0.5, 0.1) - "die meisten Münzen sind fair"
# Daten: 3 Würfe, 2x Kopf

rm(list=ls(all=T))

# 1. Definiere Parameter
prior_mean <- 0.5
prior_sd <- 0.1
n_trials <- 3
n_heads <- 2

# 2. Erstelle Gitter für mögliche p-Werte (Wahrscheinlichkeit für Kopf)
p_values <- seq(0.01, 0.99, length.out = 1000)
p_values

# 3. Berechne Prior-Dichte
prior_density <- dnorm(p_values, prior_mean, prior_sd)

# 4. Berechne Likelihood (Binomialverteilung) -> Likelihood: Binomial(2 Erfolge, 3 Versuche) für jeden möglichen p-Wert.
likelihood <- dbinom(n_heads, size = n_trials, prob = p_values)

# Die Likelihood ist NICHT die Wahrscheinlichkeit für p! Es ist die Wahrscheinlichkeit für die Daten bei gegebenem p.
# Likelihood sagt: "p = 0.67 erklärt unsere Daten (2 von 3 Kopf) am besten"
# Prior sagt: "Aber p = 0.5 ist a priori wahrscheinlicher (faire Münzen sind häufiger)"
# Posterior kombiniert beide: "Wahrscheinlich liegt p irgendwo zwischen 0.5 und 0.67"

# 5. Berechne unnormalisierte Posterior
posterior_unnorm <- prior_density * likelihood # Bayes Theorem! Posterior ∝ Prior × Likelihood 
hist(posterior_unnorm, breaks = 100)
# Problem: Die Werte in posterior_unnorm summieren sich nicht zu 1, sind also keine gültige Wahrscheinlichkeitsverteilung.

# 6. Normalisierung des Posterior / numerische Integration
posterior <- posterior_unnorm / sum(posterior_unnorm * (p_values[2] - p_values[1])) # Riemannsche Summe - eine Näherung des Integrals

# 7. Visualisierung
par(mfrow = c(2, 2)) # Vorbereitung panel space
# Prior
plot(p_values, prior_density, type = "l", col = "blue", lwd = 2,
     main = "Prior: N(0.5, 0.1)", xlab = "p (Wahrscheinlichkeit Kopf)", ylab = "Dichte")

# Likelihood
plot(p_values, likelihood, type = "l", col = "red", lwd = 2,
     main = "Likelihood: 2 Kopf aus 3 Würfen", xlab = "p", ylab = "Likelihood")

# Posterior
plot(p_values, posterior, type = "l", col = "purple", lwd = 2,
     main = "Posterior", xlab = "p", ylab = "Dichte")

# Vergleich aller drei
plot(p_values, prior_density/max(prior_density), type = "l", col = "blue", lwd = 2,
     main = "Vergleich (normalisiert)", xlab = "p", ylab = "Relative Dichte",
     ylim = c(0, 1.1))
lines(p_values, likelihood/max(likelihood), col = "red", lwd = 2)
lines(p_values, posterior/max(posterior), col = "purple", lwd = 2)


# 8. Berechne Posterior-Statistiken
posterior_mean <- sum(p_values * posterior * (p_values[2] - p_values[1]))
posterior_var <- sum((p_values - posterior_mean)^2 * posterior * (p_values[2] - p_values[1]))
posterior_sd <- sqrt(posterior_var)

# 9. Finde Maximum a Posteriori (MAP)
map_index <- which.max(posterior)
map_estimate <- p_values[map_index]

# 10. Berechne Credible Interval (95%)
cumulative_posterior <- cumsum(posterior * (p_values[2] - p_values[1]))
ci_lower <- p_values[which(cumulative_posterior >= 0.025)[1]]
ci_upper <- p_values[which(cumulative_posterior >= 0.975)[1]]

# Ergebnisse ausgeben
cat("=== POSTERIOR-STATISTIKEN ===\n")
cat("Posterior Mean:", round(posterior_mean, 3), "\n")
cat("Posterior SD:", round(posterior_sd, 3), "\n")
cat("MAP Estimate:", round(map_estimate, 3), "\n")
cat("95% Credible Interval: [", round(ci_lower, 3), ",", round(ci_upper, 3), "]\n")

# 11. Vorhersage für nächsten Wurf
cat("\n=== VORHERSAGE FÜR NÄCHSTEN WURF ===\n")
prob_next_heads <- posterior_mean
cat("Wahrscheinlichkeit für Kopf beim nächsten Wurf:", round(prob_next_heads, 3), "\n")

# Reset plot layout
par(mfrow = c(1, 1))

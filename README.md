# Quantifying uncertainty: A comparison of Bayesian and frequentist meta-analytic frameworks applied to lead exposure and children's IQ

**Authors:** Paulina Sell, Margaux Sanchez, Philippe Palmont, Simon Steiger and Dietrich Plass\
**Journal:** *[Journal name]* · [Year] · DOI: `[doi]`\
**Funding:** European Partnership for the Assessment of Risks from Chemicals (PARC) - EU Horizon (Grant Agreement No 101057014)

[![DOI](https://img.shields.io/badge/DOI-[doi]-blue)](https://doi.org/[doi])
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC--BY--4.0-lightgrey)](https://creativecommons.org/licenses/by/4.0/)
![](manuscript/tables_figures/main/graphical_abstract.png)

## About
This repository contains all code and data needed to reproduce the analyses in the paper above.
The paper compares Bayesian and frequentist meta-analytic frameworks, using lead exposure and children's IQ as a case study, in order to evaluate differences in uncertainty quantification and applicability for risk assessment

## Repository overview

- `code/` – Main R analysis scripts, including visualisations
- `data/` – Input datasets used in the meta-analysis
- `manuscript/` – Quarto files and more to reproduce the manuscript, containing rendered manuscript
- `manuscript/tables_figures/main` – All tables and figures used in main article
- `manuscript/tables_figures/supplement` – All tables and figures used in main supplementary files
- `models/` – .rds files containing model fits for main and sensitivity analyses
- `results/`– .csv files containing extracted draws from brms model fits
- `renv.lock` – R package versions for reproducibility
- `.gitignore` – Files not tracked by git
- `README.md` – This file

## Reproducing the analysis

### Requirements
- R ≥ 4.4.2
- Key packages: `metafor` (4.8.0), `brms` (2.22.0), `tidyverse` (2.0.0), `tidybayes` (3.0.7), `HDInterval`(0.2.4), `posterior` (1.6.1)
- Because `brms` is based on Stan, a C++ compiler is required. The software Rtools (available on https://cran.r-project.org/bin/windows/Rtools/) comes with a C++ compiler for Windows. On Mac, you should install Xcode. If you don't want to install additional software, you can still run most of the code, because it uses saved RDS files.
- Package versions are locked via [`renv`](https://rstudio.github.io/renv/). Run the following in R to restore the exact environment:

```r
install.packages("renv")
renv::restore()
```
> [!NOTE]
> **Environment restore may take 15–30 minutes**. If you only want to browse the code or use the saved model fits (.rds files), you can skip `renv::restore()` and install packages manually as needed.

### Running the code to reproduce the manuscript
Open `manuscript/manuscript.qmd` in RStudio or Positron and click **Render**.
The manuscript and supplementary files read pre-computed model fits (`.rds`) and figures from `models/` and `manuscript/tables_figures/`. To regenerate these from scratch instead of relying on the checked-in fits, run the scripts in `code/` in the order below — each stage depends on `.rds` files produced by the previous one.

**Stage 1 — main model fits**
1. `code/Main-model_log.R` — fits the main Bayesian model. The `saveRDS()` calls that write `models/main/m.brm_full.rds` and `models/main/fitPrior_full.rds` are commented out by default (the script re-reads the checked-in fits instead of overwriting them); uncomment them to persist a fresh fit.
2. `code/freq_model_log.R` — fits the main frequentist model (`test = "knha"` small-sample adjustment) and writes `models/main/freq_full.rds`. This is the single source for the frequentist estimate; every downstream script and `.qmd` file reads this `.rds` rather than refitting.

**Stage 2 — sensitivity model fits** (needs Stage 1's `.rds` files)
3. `code/sensitivity_analyses.R` — fits the risk-of-bias subset models (Bayesian + frequentist, both with `knha`) and the main/wide/narrow prior variants; writes `models/sensitivity/RoB/*.rds`, `models/sensitivity/priors/*.rds`, and `manuscript/tables_figures/main/table_sensi.rds`. Ignore the leave-one-out section at the bottom of this script (`results/loo_freq_results.rds`, `results/loo_bayes_results.rds`) — it is obsolete and no longer feeds the manuscript.
4. `code/sensitivity_priors.R` (sources `code/prior_configs.R`) — fits the one-way prior sensitivity models; some configs reuse Stage 1/3 fits (`m.brm_full.rds`, `priors/{narrow,wide}`) rather than refitting. Writes `models/sensitivity/priors_oneway/*.rds` and `manuscript/tables_figures/main/table_sensi_oneway.rds`.

**Stage 3 — figures** (need Stage 1/2 `.rds` files)
5. `code/Bayes_forest_plot.R` — combined Bayesian + frequentist forest plot, `manuscript/tables_figures/main/figure2.png`.
6. `code/visualisations.R` — remaining main and supplementary figures.
7. `code/figures_priors_oneway.R` (sources `code/prior_configs.R`) — prior-comparison figures for the supplement.
8. `code/transformations_vis.R` is a standalone exploratory plot (log vs. linear exposure-response), not wired into any manuscript output.

> [!NOTE]
> `results/traceplot.png`, embedded in `file5_additional-plots.qmd`, currently has no generating script in `code/` — it must be supplied manually.

**Stage 4 — render Quarto files** (need Stage 1–3 outputs)
- `manuscript/manuscript.qmd` — main manuscript. Also embeds `tables_figures/main/figure1.png` (PRISMA diagram), which is created manually, not by code.
- `manuscript/tables_figures/supplement/model_output.qmd`
- `manuscript/tables_figures/supplement/file5_additional-plots.qmd`
- `manuscript/tables_figures/supplement/I2_calculation.qmd` — text only, no code dependencies, can be rendered independently at any time.

`manuscript/tables_figures/supplement/loo_analysis.qmd` is obsolete and no longer part of the manuscript; it is not included in the render order above.

## Resources

Resources I found helpful while working on this project:

- [Harrer, M. et al. (2021). *Doing Meta-Analysis With R: A Hands-On Guide.*](https://doing-meta.guide/bayesian-ma) — accessible introduction to Bayesian meta-analysis in R
- [Grant, R., & Di Tanna, G. L. (2025). *Bayesian Meta-Analysis: A Practical Introduction.*](https://doi.org/10.1201/9781003375821) - great book on Bayesian meta-analysis

## Citation

```bibtex
@article{sell2026,
  author  = {[Last, First] and [Last, First] and …},
  title   = {[Paper title]},
  journal = {[Journal name]},
  year    = {[Year]},
  doi     = {[doi]}
}
```

## License

This repository is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). 
You are free to share and adapt the material for any purpose, provided appropriate credit is given.

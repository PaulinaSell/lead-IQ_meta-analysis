# Quantifying uncertainty: A comparison of Bayesian and frequentist meta-analytic frameworks applied to lead exposure and children's IQ

**Authors:** Paulina Sell, Margaux Sanchez, Philippe Palmont, Simon Steiger and Dietrich Plass
**Journal:** *[Journal name]* · [Year] · DOI: `[doi]`  
**Funding:** European Partnership for the Assessment of Risks from Chemicals (PARC) - EU Horizon (Grant Agreement No 101057014)

[![DOI](https://img.shields.io/badge/DOI-[doi]-blue)](https://doi.org/[doi])
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC--BY--4.0-lightgrey)](LICENSE)

![](manuscript/tables_figures/main/graphical_abstract.png){width=100% height=100% fig-align="center"}

This repository contains all code and data needed to reproduce the analyses in the paper above.
The paper compares Bayesian and frequentist meta-analytic frameworks, using lead exposure and children's IQ as a case study, in order to evaluate differences in uncertainty quantification and applicability for risk assessment

| Path | Description |
|---|---|
| `code/` | Main R analysis scripts, including visualisations |
| `data/` | Input datasets used in the meta-analysis |
| `manuscript/` | Quarto files and more to reproduce the manuscript, containing rendered manuscript |
| `manuscript/tables_figures/main` | All tables and figures used in main article |
| `manuscript/tables_figures/supplement` | All tables and figures used in main supplementary files |
| `models/` | RDS files containing model fits for main and sensitivity analyses |
| `results/`| CSV files containing extracted draws from brms model fits |
| `renv.lock` | R package versions for reproducibility |
| `.gitignore` | Files not tracked by git |
| `README.md` | This file |

## Reproducing the analysis

### Requirements
- R (≥ 4.4.2 (2024-10-31))
- Key packages: `metafor` (4.8.0), `brms` (2.22.0), `tidyverse` (2.0.0), `tidybayes` (3.0.7), `HDInterval`(0.2.4), `posterior` (1.6.1)
- Package versions are locked via [`renv`](https://rstudio.github.io/renv/). Run the following in R to restore the exact environment:

```r
install.packages("renv")
renv::restore()
```

### Running the code to reproduce the manuscript
Open `manuscript/manuscript.qmd` in RStudio and click **Render**

## Resources

Resources I found helpful while working on this project:

- [Harrer, M. et al. (2021). Doing Meta-Analysis With R: A Hands-On Guide.](https://doing-meta.guide/bayesian-ma) — accessible introduction to Bayesian meta-analysis in R



- [Viechtbauer (2010)](https://doi.org/[doi]) — the paper introducing the `metafor` package

**Bayesian methods in R**
- [Bürkner (2017)](https://doi.org/[doi]) — the paper introducing the `brms` package
- [McElreath (2020) *Statistical Rethinking*](https://doi.org/[doi]) — excellent conceptual introduction to Bayesian statistics

**Lead exposure and children's IQ**
- [Source](url) — [one-line description]
- [Source](url) — [one-line description]

## Citation

```bibtex
@article{sell[year],
  author  = {[Last, First] and [Last, First] and …},
  title   = {[Paper title]},
  journal = {[Journal name]},
  year    = {[Year]},
  doi     = {[doi]}
}
```

## License

This repository is licensed under [CC BY 4.0](LICENSE). You are free to share and adapt the material for any purpose, provided appropriate credit is given.

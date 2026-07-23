# Quick GLMM Review

An interactive R Shiny app for fitting and checking **generalized linear mixed models** with [`glmmTMB`](https://cran.r-project.org/package=glmmTMB) — a point-and-click companion for exploratory GLMM work, with proper simulation-based residual diagnostics via [`DHARMa`](https://cran.r-project.org/package=DHARMa) and post-hoc comparisons via [`emmeans`](https://cran.r-project.org/package=emmeans).

It's a sibling to a companion `lmerTest`-based app for linear mixed models, extended to counts, proportions, zero-inflated data, and binary outcomes.

## Features

- **Two independent tabs**
  - **General GLMM** — continuous, count, and proportion responses
  - **Binary (0/1) GLMM** — kept separate, with a restricted family list, link-function picker (logit / probit / cloglog / cauchit), and no dispersion model, since Bernoulli data doesn't have a free dispersion parameter
- **Family picker**: Gaussian, Gamma, Poisson, negative binomial (nbinom1 / nbinom2), Tweedie, Beta
- **Zero-inflation (`ziformula`)**: toggle on/off, optionally add predictors — builds `~0`, `~1`, or `~var1 + var2`
- **Dispersion model (`dispformula`)**: 0, 1, or 2 predictors, including any user-created interaction/combined variable
- **Random effects**: multiple grouping factors, optional random slope attached to a chosen group
- **Diagnostics**: DHARMa simulated-residual QQ/residual plot, overdispersion test, zero-inflation test, outlier test
- **Inference**: Type III Wald chi-square tests (`car::Anova`), EMMeans, pairwise comparisons, compact letter display, and an EMMeans plot
- **Data tools**: CSV upload or built-in example dataset, numeric-to-factor overrides, and a combined/interaction-variable creator
- **Reproducibility**: downloadable R code, CSVs, and plots for everything the app produces

## Requirements

```r
install.packages(c(
  "shiny", "shinyjs", "bslib", "glmmTMB", "DHARMa", "emmeans",
  "multcomp", "multcompView", "ggplot2", "DT", "car"
))
```

R ≥ 4.1 is recommended (matches `glmmTMB`'s current requirements).

## Running the app

```r
shiny::runApp("glmmTMB_app.R")
```

or open `glmmTMB_app.R` in RStudio and click **Run App**.

## Branding

The app looks for an optional `components.R` file in the same directory (shared UF/IFAS theme: colors, logo, tooltip helper, clipboard-copy JS). If it isn't present, the app falls back to an inline theme and still runs standalone — `components.R` is not required.

## Using your own data

1. Choose **Upload a CSV file** as the data source.
2. Use **Variable types** to mark any numeric ID/treatment columns (e.g. `Block`, `Nitrogen`) as factors if their *levels* should be compared, not fit as a continuous slope.
3. Optionally create a combined/interaction variable (e.g. `Treatment.Season`) if you want to use it as a dispersion-model predictor.
4. Pick a response, family, fixed/random effects, and (optionally) zero-inflation and dispersion terms, then **Fit model**.
5. Check the **DHARMa residuals** tab before trusting the **Wald ANOVA** or **EMMeans** results.

## A note on binary outcomes

0/1 data is intentionally handled on its own tab rather than folded into the general family dropdown. Binomial/Bernoulli responses don't have a meaningfully estimable dispersion parameter the way counts or continuous data do, and DHARMa's residual plot is naturally grainy for binary data — the formal tests matter more than the visual pattern there. Keeping the tab separate avoids silently misconfiguring a general model for a binary response.

## Known limitations

- No automatic detection of which family best fits the data — try a couple and compare AIC/DHARMa diagnostics.
- Kenward-Roger/Satterthwaite-style denominator df (available for `lmerTest` models) aren't used here; inference uses Wald chi-square tests instead, which is the standard approach for `glmmTMB`.
- Very small or unbalanced datasets may cause `glmmTMB` convergence warnings — these surface in the **Model & code → Model summary** panel.

## License

Add your preferred license here (e.g. MIT).

# crown: Cluster-Randomized Trial Analysis with Outcomes Weighted for Nonresponse

`crown` combines a cluster-randomized trial with auxiliary baseline data to
estimate population mean potential outcomes and population-average causal
effects on the risk-difference and risk-ratio scales.

## Installation

From the folder containing `DESCRIPTION`:

```r
install.packages(c("xgboost", "SuperLearner", "numDeriv"))
install.packages(".", repos = NULL, type = "source")
library(crown)
```

The package requires R >= 4.1 and XGBoost >= 3.2.0. The source code is in the
[PopART repository](https://github.com/zhitanhe21/PopART).

## Quick start

The package includes a weighted example with 50,000 trial observations and
50,000 auxiliary observations:

```r
library(crown)
data(crown_example)

trial <- crown_example[crown_example$S == 1, ]
auxiliary <- crown_example[crown_example$S == 0, ]
```

Run the default analysis, Proposed AIPW with cross-fitted XGBoost:

```r
fit <- crown(
  trial_data = trial,
  auxiliary_data = auxiliary,
  outcome = "Y",
  treatment = "A",
  response = "R",
  censoring = "C",
  cluster = "cluster",
  covariates = c("X1", "X2", "W1", "W2"),
  version = "proposed",
  estimator = "aipw",
  model = "xgboost",
  auxiliary_weight = "sampling_weight"
)

print(fit)
summary(fit)
```

The complete script is also installed with the package:

```r
source(system.file("examples", "run_analysis.R", package = "crown"))
```

### Example result

| Parameter | Estimate | Standard error | 95% CI |
|---|---:|---:|---|
| `eta(0)` | 0.285 | 0.005 | [0.275, 0.294] |
| `eta(1)` | 0.427 | 0.006 | [0.415, 0.438] |
| `RD` | 0.142 | 0.007 | [0.127, 0.156] |
| `RR` | 1.498 | 0.032 | [1.436, 1.560] |

This run took 5.44 seconds with R 4.5.2 and XGBoost 3.2.1.1.
These are results from one synthetic dataset, not Monte Carlo averages.

## Choose the analysis

The main choices are made directly in `crown()`:

| Argument | Default | Choices |
|---|---|---|
| `version` | `"proposed"` | `"proposed"`, `"naive"` |
| `estimator` | `"aipw"` | `"aipw"`, `"gformula"`, `"ipw"` |
| `model` | `"xgboost"` | `"xgboost"`, `"logistic"` |

`proposed` uses the auxiliary sample to represent the target population.
`naive` uses the trial responders and does not adjust for nonresponse.

All 12 combinations return point estimates. Standard errors and 95% confidence
intervals are currently available for:

| Model | Version | Estimator | CI |
|---|---|---|---|
| Logistic | Naive or Proposed | G-formula, IPW, or AIPW | Yes |
| XGBoost | Proposed | AIPW | Yes |
| XGBoost | Proposed | G-formula or IPW | No |
| XGBoost | Naive | G-formula, IPW, or AIPW | No |

XGBoost uses cross-fitting. `K = 5L` sets the number of folds. Use `arguments`
to pass XGBoost settings and `random_seed` to reproduce the fold split.

## Prepare your data

Pass the trial and auxiliary samples as separate data frames. One row represents
one person.

| Variable | Trial | Auxiliary |
|---|---|---|
| Cluster ID | Required | Required; must match a trial cluster |
| Treatment | Required and constant within cluster | Inferred from the trial |
| Response | Required | Not required |
| Censoring | Required | Not required |
| Outcome | Required for uncensored responders | Not required |
| Covariates | Required | The same covariates are required |
| Sampling weight | Not used | Optional |

Use the column names when calling `crown()`. Missing outcomes are allowed for
nonresponders and censored responders.

### Auxiliary sampling weights

If the auxiliary sample already represents the target population, leave the
weight unspecified:

```r
auxiliary_weight = NULL
```

If it has a sampling-weight column, provide its name:

```r
auxiliary_weight = "sampling_weight"
```

The weights are normalized to have mean one. They are used only by Proposed
analyses; Naive analyses do not use the auxiliary sample.

### Included weighted data

`crown_example` has 100,000 rows and 200 clusters. Each cluster contains 250
trial and 250 auxiliary observations. The auxiliary sample contains 150 people
with `W1 = 1` and 100 with `W1 = 0` in each cluster, so `W1 = 1` is deliberately
oversampled. `sampling_weight` corrects this unequal sampling. The compressed
dataset is about 364 KB.

## Main functions

| Function | What it does |
|---|---|
| `crown()` | Runs one selected analysis |
| `fit_gformula()` | Fits parametric G-formula |
| `fit_ipw()` | Fits parametric IPW |
| `fit_aipw()` | Fits parametric AIPW |
| `dml_fit()` | Fits the lower-level cross-fitted estimator |

Use `print(fit)` for the four main estimates and `summary(fit)` for estimates,
standard errors, and confidence intervals. Unavailable intervals are shown as
`NA`.

## Simulations

The simulation code stays in the GitHub repository and is not included in the
installed package. From the `crown` source folder, run:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "ggh4x", "legendry"))
source("simulation/sim_scripts/ss1.R")
source("simulation/sim_scripts/ss2.R")
```

- Simulation 1 compares Naive and Proposed parametric G-formula, IPW, and AIPW.
- Simulation 2 compares Proposed parametric AIPW with cross-fitted XGBoost AIPW.
- Both use `(500, 500)`, `(500, 5000)`, `(5000, 500)`, and `(5000, 5000)` for
  `(n_trial, n_auxiliary)`.
- Both currently run 20 Monte Carlo replicates per sample-size combination.
- Replicates and folds run sequentially; XGBoost may use its own threads.

The figures below come from earlier 10,000-replicate runs and still need to be
verified against the current code.

### Simulation 1

![Simulation 1 risk-difference estimates](simulation/sim_figures/reference/sim1_estimates.png)

![Simulation 1 variance comparison](simulation/sim_figures/reference/sim1_variance.png)

![Simulation 1 confidence-interval coverage](simulation/sim_figures/reference/sim1_confidence.png)

### Simulation 2

![Simulation 2 risk-difference estimates](simulation/sim_figures/reference/sim2_estimates.png)

![Simulation 2 variance comparison](simulation/sim_figures/reference/sim2_variance.png)

![Simulation 2 confidence-interval coverage](simulation/sim_figures/reference/sim2_confidence.png)

## References

To be added.

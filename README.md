# crown

`crown` combines a cluster-randomized trial with an auxiliary sample to
estimate population mean potential outcomes, risk differences, and risk
ratios when trial outcomes are affected by nonresponse or censoring.

## Installation

Install the development version from GitHub:

```r
install.packages("devtools")
library(devtools)
install_github("brian-d-richardson/crown", ref = "crown-1.0.0")
library(crown)
```

## Supported analyses

| Analysis | `version` | `estimator` | `model` |
|---|---|---|---|
| Naive G-formula | `"naive"` | `"gformula"` | `"logistic"` |
| Proposed G-formula | `"proposed"` | `"gformula"` | `"logistic"` |
| Naive IPW | `"naive"` | `"ipw"` | `"logistic"` |
| Proposed IPW | `"proposed"` | `"ipw"` | `"logistic"` |
| Naive AIPW | `"naive"` | `"aipw"` | `"logistic"` |
| Proposed AIPW | `"proposed"` | `"aipw"` | `"logistic"` |
| Proposed DML | `"proposed"` | `"aipw"` | `"xgboost"` |

These are the only supported combinations.

## Data requirements

| Variable | Trial | Auxiliary |
|---|---|---|
| Cluster ID | Required | Required; must match a trial cluster |
| Treatment | Required and constant within cluster | Not required |
| Response indicator | Required | Not required |
| Censoring indicator | Required | Not required |
| Outcome | Required for uncensored responders | Not required |
| Covariates | Required | Same covariates required |
| Sampling weight | Not used | Optional |

Pass the trial and auxiliary samples as separate data frames. If auxiliary
sampling weights are unequal, provide the weight-column name through
`auxiliary_weight`; otherwise leave it as `NULL`.

## Application

The package includes `crown_example`, a synthetic dataset with 50,000 trial
participants and 50,000 auxiliary participants across 20 clusters. Its first
five participants are:

| id | cluster | S | A | R | C | Y | X1 | X2 | W1 | W2 | sampling_weight |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 11 | 1 | 1 | 1 | 0 | 1 | -2.1958 | 1 | 1 | 0.3342 | 1 |
| 2 | 11 | 1 | 1 | 1 | 0 | 0 | -2.1958 | 1 | 0 | 0.8751 | 1 |
| 3 | 11 | 1 | 1 | 1 | 0 | 0 | -2.1958 | 1 | 1 | -0.0489 | 1 |
| 4 | 11 | 1 | 1 | 0 | NA | NA | -2.1958 | 1 | 1 | -0.9554 | 1 |
| 5 | 11 | 1 | 1 | 1 | 0 | 0 | -2.1958 | 1 | 0 | -0.2128 | 1 |

Select an analysis through `version`, `estimator`, and `model`:

```r
data(crown_example)

trial <- crown_example[crown_example$S == 1, ]
auxiliary <- crown_example[crown_example$S == 0, ]

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

summary(fit)
```

This call fits Proposed DML and returns:

| Parameter | Estimate | Standard error | 95% CI |
|---|---:|---:|---:|
| Mean under control | 0.272 | 0.005 | [0.263, 0.281] |
| Mean under treatment | 0.404 | 0.006 | [0.392, 0.416] |
| Risk difference | 0.132 | 0.007 | [0.118, 0.147] |
| Risk ratio | 1.486 | 0.033 | [1.422, 1.550] |

## Simulations

The simulation scripts run 50 Monte Carlo replicates per sample-size setting by
default. Change `mc_reps` near the top of each script for a different run:

```r
source("simulation/sim_scripts/ss1.R")
source("simulation/sim_scripts/ss2.R")
```

Each script saves `sd_local.csv` and produces point-estimate, variance, and
confidence-interval coverage figures. The figures below use the 10,000
replicates per setting stored in `simulation/sim_data/`.

### Simulation 1

Simulation 1 uses a 20-cluster randomized trial and an auxiliary sample with a
different covariate distribution. Response, censoring, and binary outcomes
depend on treatment and covariates. Correct and misspecified outcome and
selection models are compared across trial and auxiliary sample sizes of 500
and 5,000.

<p align="center"><strong>Parameter estimates</strong><br>
  <img src="simulation/sim_figures/sim1/sim1_estimates.png"
       alt="Simulation 1 estimates" width="760">
</p>

<p align="center"><strong>Variance estimates</strong><br>
  <img src="simulation/sim_figures/sim1/sim1_variance.png"
       alt="Simulation 1 variances" width="760">
</p>

<p align="center"><strong>Confidence-interval coverage</strong><br>
  <img src="simulation/sim_figures/sim1/sim1_confidence.png"
       alt="Simulation 1 confidence-interval coverage" width="760">
</p>

### Simulation 2

Simulation 2 uses the same cluster and sample-size settings but makes the
response and outcome mechanisms nonlinear through squared and sinusoidal
covariate effects. It evaluates all seven supported estimators; the figures
show the four Proposed estimators.

<p align="center"><strong>Parameter estimates</strong><br>
  <img src="simulation/sim_figures/sim2/sim2_estimates.png"
       alt="Simulation 2 estimates" width="1000">
</p>

<p align="center"><strong>Variance estimates</strong><br>
  <img src="simulation/sim_figures/sim2/sim2_variance.png"
       alt="Simulation 2 variances" width="1000">
</p>

<p align="center"><strong>Confidence-interval coverage</strong><br>
  <img src="simulation/sim_figures/sim2/sim2_confidence.png"
       alt="Simulation 2 confidence-interval coverage" width="1000">
</p>

## Reference

To be added after publication.

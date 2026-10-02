# Import data

library(devtools)

package_dir <- normalizePath(
  file.path(dirname(sys.frame(1)$ofile), "..", "..")
)
load_all(package_dir)

data(crown_example)
print(head(crown_example, 5), row.names = FALSE)

# Split trial and auxiliary samples
# The data contain 20 clusters and 100,000 observations in total.

trial <- crown_example[crown_example$S == 1, ]
auxiliary <- crown_example[crown_example$S == 0, ]

# Fit Crown
# Supported combinations:
# Logistic: naive/proposed G-formula, IPW, and AIPW.
# XGBoost: proposed AIPW (DML) only.

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

# Results

print(fit)
print(summary(fit))

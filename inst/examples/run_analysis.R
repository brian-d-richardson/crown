# Import data

library(crown)

data(crown_example)

# Split trial and auxiliary samples

trial <- crown_example[crown_example$S == 1, ]
auxiliary <- crown_example[crown_example$S == 0, ]

# Fit Crown

fit <- crown(
  trial_data = trial,
  auxiliary_data = auxiliary,
  outcome = "Y",
  treatment = "A",
  response = "R",
  censoring = "C",
  cluster = "cluster",
  covariates = c("X1", "X2", "W1", "W2"),
  version = "proposed",       # "proposed" or "naive"
  estimator = "aipw",        # "aipw", "gformula", or "ipw"
  model = "xgboost",         # "xgboost" or "logistic"
  auxiliary_weight = "sampling_weight"
)

# Results

print(fit)
print(summary(fit))

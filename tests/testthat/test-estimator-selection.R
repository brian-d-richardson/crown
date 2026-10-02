test_that("the seven supported combinations return complete inference", {
  d <- example_data()
  labels <- c(aipw = "AIPW", gformula = "G-Formula", ipw = "IPW")
  supported <- rbind(
    expand.grid(
      version = c("naive", "proposed"),
      estimator = names(labels),
      model = "logistic",
      stringsAsFactors = FALSE
    ),
    data.frame(version = "proposed", estimator = "aipw", model = "xgboost")
  )

  for (i in seq_len(nrow(supported))) {
    setting <- supported[i, ]
    expect_warning(fit <- crown(subset(d, S == 1), subset(d, S == 0),
      "Y", "A", "R", "C", "cluster", c("X", "W"),
      version = setting$version, model = setting$model,
      estimator = setting$estimator, K = 2L,
      arguments = list(nrounds = 3L)), NA)
    expect_equal(nrow(fit$estimates), 4L)
    expect_setequal(fit$estimates$estimator, labels[[setting$estimator]])
    expect_named(fit$covariance,
                 paste(setting$estimator, setting$version, sep = "_"))
    expect_true(all(is.finite(fit$estimates$estimate)))
    expect_true(all(is.finite(fit$estimates$std_error)))
  }
})

test_that("unsupported XGBoost combinations are rejected", {
  d <- example_data()
  unsupported <- expand.grid(
    version = c("naive", "proposed"),
    estimator = c("aipw", "gformula", "ipw"),
    stringsAsFactors = FALSE
  )
  unsupported <- unsupported[
    unsupported$version != "proposed" | unsupported$estimator != "aipw", ]

  for (i in seq_len(nrow(unsupported))) {
    setting <- unsupported[i, ]
    expect_error(crown(subset(d, S == 1), subset(d, S == 0),
      "Y", "A", "R", "C", "cluster", c("X", "W"),
      version = setting$version, model = "xgboost",
      estimator = setting$estimator),
      "XGBoost is supported only for Proposed AIPW")
  }
})

test_that("DML fits four nuisance models in each fold", {
  calls <- 0L
  local_mocked_bindings(
    nonpar_est = function(x, y, arguments, wts = NULL) {
      calls <<- calls + 1L
      NULL
    },
    nonpar_pred = function(mod, newdata) rep(.5, nrow(newdata))
  )
  d <- example_data()
  crown(subset(d, S == 1), subset(d, S == 0),
    "Y", "A", "R", "C", "cluster", c("X", "W"), K = 2L)
  expect_equal(calls, 8L)
})

test_that("logistic runs only the requested fitting function", {
  calls <- character()
  row <- .eta_result(.2, .3, diag(.01, 2L))
  local_mocked_bindings(
    fit_aipw = function(...) { calls <<- c(calls, "aipw"); row },
    fit_gformula = function(...) { calls <<- c(calls, "gformula"); row },
    fit_ipw = function(...) { calls <<- c(calls, "ipw"); row }
  )
  d <- example_data()
  for (estimator in c("aipw", "gformula", "ipw")) {
    calls <- character()
    crown(subset(d, S == 1), subset(d, S == 0),
      "Y", "A", "R", "C", "cluster", c("X", "W"),
      model = "logistic", estimator = estimator)
    expect_identical(calls, estimator)
  }
})

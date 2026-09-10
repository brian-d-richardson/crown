test_that("the logistic interface defaults to proposed AIPW only", {
  d <- example_data()
  local_mocked_bindings(dml_fit = function(...) stop("DML should not run"))
  fit <- crown(subset(d, S == 1), subset(d, S == 0),
    "Y", "A", "R", "C", "cluster", c("X", "W"), model = "logistic")
  expect_s3_class(fit, "crown_fit")
  expect_equal(nrow(fit$estimates), 4L)
  expect_length(fit$covariance, 1L)
  expect_setequal(fit$estimates$estimator, "AIPW")
  expect_setequal(fit$estimates$version, "Proposed")
  expect_null(fit$dml)
  expect_true(all(is.finite(fit$estimates$std_error)))
  expect_output(print(fit), "Crown causal estimates")
  expect_false(any(grepl("cross-fitting", capture.output(print(fit)), fixed = TRUE)))
  expect_equal(nrow(summary(fit)), 4L)
  for (name in names(fit$covariance)) {
    expect_equal(fit$covariance[[name]], t(fit$covariance[[name]]))
  }
})

test_that("each logistic version agrees with direct fits and uses its formulas", {
  d <- example_data()
  d$wt[d$S == 0] <- seq(.5, 1.5, length.out = sum(d$S == 0))
  for (version in c("naive", "proposed")) {
    weight_formula <- if (version == "naive") C ~ W else Q ~ W
    rows <- list(
      fit_gformula(d, Y ~ A * W, version),
      fit_ipw(d, weight_formula, version),
      fit_aipw(d, Y ~ A * W, weight_formula, version)
    )
    names(rows) <- paste(c("gformula", "ipw", "aipw"), version, sep = "_")
    for (estimator in c("gformula", "ipw", "aipw")) {
      fit <- crown(subset(d, S == 1), subset(d, S == 0),
        "Y", "A", "R", "C", "cluster", c("X", "W"),
        version = version, model = "logistic", auxiliary_weight = "wt",
        outcome_formula = Y ~ A * W, propensity_formula = Q ~ W,
        censoring_formula = C ~ W, estimator = estimator)
      expected <- .format_crown_results(rows[paste(estimator, version, sep = "_")])
      expect_equal(fit$estimates, expected$estimates)
      expect_equal(fit$covariance, expected$covariance)
    }
  }
})

test_that("the default interface runs only proposed cross-fitted XGBoost DML", {
  d <- example_data()
  d[d$S == 1 & d$R == 0, c("X", "W")] <- NA_real_
  combined <- .prepare_crown_data(subset(d, S == 1), subset(d, S == 0),
    "Y", "A", "R", "C", "cluster", c("X", "W"), NULL, c(0, 1))
  expected <- dml_fit(combined, c("X", "W"), c("X", "W"), K = 2L,
    arguments = list(nrounds = 3L), random_seed = 10L)
  local_mocked_bindings(
    fit_gformula = function(...) stop("G-formula should not run"),
    fit_ipw = function(...) stop("IPW should not run"),
    fit_aipw = function(...) stop("Parametric AIPW should not run")
  )
  fit <- crown(subset(d, S == 1), subset(d, S == 0),
    "Y", "A", "R", "C", "cluster", c("X", "W"),
    K = 2L,
    arguments = list(nrounds = 3L), random_seed = 10L)
  explicit <- crown(subset(d, S == 1), subset(d, S == 0),
    "Y", "A", "R", "C", "cluster", c("X", "W"),
    version = "proposed", model = "xgboost", estimator = "aipw", K = 2L,
    arguments = list(nrounds = 3L), random_seed = 10L)
  expect_equal(fit, explicit)
  expect_equal(fit$dml, expected)
  expect_equal(nrow(fit$estimates), 4L)
  expect_length(fit$covariance, 1L)
  expect_setequal(fit$estimates$estimator, "AIPW")
  expect_setequal(fit$estimates$version, "Proposed")
  expect_output(print(fit), "AIPW")
  expect_output(print(fit), "2 folds of cross-fitting used", fixed = TRUE)
  expect_false(any(grepl("DML", capture.output(print(fit)), fixed = TRUE)))
  expect_setequal(summary(fit)$estimator, "AIPW")
  expect_named(fit$covariance, "aipw_proposed")
  expected_row <- .eta_result(expected$eta_hat[1], expected$eta_hat[2],
                             expected$eta_hat_cov)
  expect_equal(fit$estimates,
               .format_crown_results(list(aipw_proposed = expected_row))$estimates)
  expect_equal(nrow(summary(fit)), 4L)
})

test_that("unsupported choices do not silently select a different method", {
  expect_error(crown(version = "unknown"), "arg")
  expect_error(crown(model = "unknown"), "arg")
  expect_error(crown(model = "parametric"), "arg")
  expect_error(crown(model = "nonparametric"), "arg")
  expect_error(crown(estimator = "unknown"), "arg")
})

test_that("cluster treatment is checked and carried into the auxiliary sample", {
  d <- example_data()
  trial <- subset(d, S == 1)
  auxiliary <- subset(d, S == 0)
  auxiliary$A <- NULL
  combined <- .prepare_crown_data(trial, auxiliary, "Y", "A", "R", "C",
                                  "cluster", c("X", "W"), NULL, c(0, 1))
  expect_equal(combined$A, d$A)
  trial$A[1] <- 1 - trial$A[1]
  expect_error(.prepare_crown_data(trial, auxiliary, "Y", "A", "R", "C",
    "cluster", c("X", "W"), NULL, c(0, 1)), "one treatment")
})

test_that("RD and RR variances use the risk covariance and delta method", {
  row <- .eta_result(.4, .2, matrix(c(.01, .002, .002, .02), 2))
  fit <- .format_crown_results(list(dml_proposed = row))
  expect_equal(fit$estimates$estimate, c(.4, .2, -.2, .5))
  expect_equal(fit$variance$variance[3], .026)
  expect_equal(fit$variance$variance[4], .128125)
})

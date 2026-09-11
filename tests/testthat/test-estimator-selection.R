test_that("all supported combinations return estimates and confidence intervals", {
  d <- example_data()
  labels <- c(aipw = "AIPW", gformula = "G-Formula", ipw = "IPW")

  for (version in c("naive", "proposed")) {
    for (estimator in names(labels)) {
      fit <- crown(subset(d, S == 1), subset(d, S == 0),
        "Y", "A", "R", "C", "cluster", c("X", "W"),
        version = version, model = "parametric", estimator = estimator)
      expect_equal(nrow(fit$estimates), 4L)
      expect_setequal(fit$estimates$estimator, labels[[estimator]])
      expect_named(fit$covariance, paste(estimator, version, sep = "_"))
      expect_true(all(is.finite(fit$estimates$estimate)))
      expect_true(all(is.finite(fit$estimates$std_error)))
    }
  }

  fit <- crown(subset(d, S == 1), subset(d, S == 0),
    "Y", "A", "R", "C", "cluster", c("X", "W"),
    version = "proposed", model = "nonparametric", estimator = "aipw",
    K = 2L, arguments = list(nrounds = 3L))
  expect_true(all(is.finite(fit$estimates$estimate)))
  expect_true(all(is.finite(fit$estimates$std_error)))
})

test_that("unsupported nonparametric combinations give a useful error", {
  d <- example_data()
  run <- function(version, estimator) {
    crown(subset(d, S == 1), subset(d, S == 0),
      "Y", "A", "R", "C", "cluster", c("X", "W"),
      version = version, model = "nonparametric", estimator = estimator)
  }

  unsupported <- list(
    c("proposed", "gformula"),
    c("proposed", "ipw"),
    c("naive", "gformula"),
    c("naive", "ipw"),
    c("naive", "aipw")
  )
  for (choice in unsupported) {
    expect_error(
      run(choice[[1]], choice[[2]]),
      "Unsupported combination"
    )
  }
  expect_error(run("naive", "aipw"), "model = 'parametric'")
  expect_error(run("proposed", "ipw"), "estimator = 'aipw'")
})

test_that("nonparametric AIPW fits all four nuisance models in each fold", {
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
    "Y", "A", "R", "C", "cluster", c("X", "W"),
    model = "nonparametric", K = 2L)
  expect_equal(calls, 8L)
})

test_that("parametric analyses run only the requested fitting function", {
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
      model = "parametric", estimator = estimator)
    expect_identical(calls, estimator)
  }
})

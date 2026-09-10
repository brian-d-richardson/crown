test_that("all version, model and estimator combinations select one estimator", {
  d <- example_data()
  labels <- c(aipw = "AIPW", gformula = "G-Formula", ipw = "IPW")
  for (version in c("naive", "proposed")) {
    for (model in c("logistic", "xgboost")) {
      for (estimator in names(labels)) {
        point_only <- model == "xgboost" &&
          (version == "naive" || estimator != "aipw")
        run <- function() crown(subset(d, S == 1), subset(d, S == 0),
          "Y", "A", "R", "C", "cluster", c("X", "W"),
          version = version, model = model, estimator = estimator,
          K = 2L, arguments = list(nrounds = 3L))
        if (point_only) {
          expect_warning(fit <- run(), "point estimates only")
          expect_true(all(is.na(fit$estimates$std_error)))
          expect_true(all(is.na(summary(fit)[["95% CI"]])))
          expect_output(print(fit), "Point estimates only")
        } else {
          expect_warning(fit <- run(), NA)
          expect_true(all(is.finite(fit$estimates$std_error)))
        }
        expect_equal(nrow(fit$estimates), 4L)
        expect_setequal(fit$estimates$estimator, labels[[estimator]])
        expect_named(fit$covariance, paste(estimator, version, sep = "_"))
        expect_true(all(is.finite(fit$estimates$estimate)))
      }
    }
  }
})

test_that("XGBoost G-formula and IPW match their own fold formulas", {
  d <- example_data()
  d$wt[d$S == 0] <- seq(.5, 1.5, length.out = sum(d$S == 0))
  for (version in c("naive", "proposed")) {
    for (estimator in c("gformula", "ipw")) {
      if (version == "proposed") {
        expect_warning(fit <- dml_fit(d, c("X", "W"), c("X", "W"),
          K = 2L, arguments = list(nrounds = 3L), estimator = estimator),
          "point estimates only")
      } else {
        expect_warning(fit <- .naive_xgboost(d, c("X", "W"), 2L,
          list(nrounds = 3L), 1L, estimator), "point estimates only")
      }
      risks <- matrix(0, 2L, 2L)
      for (k in 1:2) {
        x <- fit$dat[fit$dat$fold == k, ]
        for (a in 0:1) {
          if (estimator == "gformula") {
            rows <- if (version == "proposed") x$S == 0 else rep(TRUE, nrow(x))
            mu <- x[[paste0("muhat_", a)]][rows]
            weights <- if (version == "proposed") x$wt[rows] else rep(1, sum(rows))
            risks[k, a + 1L] <- sum(mu * weights) / sum(weights)
          } else {
            rows <- x$S == 1 & x$R == 1 & x$C == 0 & x$A == a
            risks[k, a + 1L] <- sum(x$Y[rows] / x$pihat[rows]) /
              sum(1 / x$pihat[rows])
          }
        }
      }
      expect_equal(unname(fit$eta_hat), colMeans(risks))
      expect_true(all(is.na(fit$eta_hat_cov)))
    }
  }
})

test_that("XGBoost skips nuisance models not needed by the estimator", {
  calls <- 0L
  local_mocked_bindings(
    nonpar_est = function(x, y, arguments, wts = NULL) {
      calls <<- calls + 1L
      NULL
    },
    nonpar_pred = function(mod, newdata) rep(.5, nrow(newdata))
  )
  d <- example_data()
  for (version in c("naive", "proposed")) {
    expected <- if (version == "naive") c(aipw = 3L, gformula = 2L, ipw = 1L) else
      c(aipw = 4L, gformula = 2L, ipw = 2L)
    for (estimator in names(expected)) {
      calls <- 0L
      suppressWarnings(crown(subset(d, S == 1), subset(d, S == 0),
        "Y", "A", "R", "C", "cluster", c("X", "W"),
        version = version, estimator = estimator, K = 2L))
      expect_equal(calls, 2L * expected[[estimator]])
    }
  }
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

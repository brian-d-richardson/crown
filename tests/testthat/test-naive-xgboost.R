test_that("naive XGBoost uses only responders and reports point estimates", {
  d <- example_data()
  run <- function(d) crown(subset(d, S == 1), subset(d, S == 0),
    "Y", "A", "R", "C", "cluster", c("X", "W"),
    version = "naive", model = "xgboost", K = 2L,
    arguments = list(nrounds = 3L), random_seed = 10L)
  expect_warning(fit <- run(d), "point estimates only")
  expect_setequal(fit$estimates$estimator, "AIPW")
  expect_setequal(fit$estimates$version, "Naive")
  expect_true(all(is.finite(fit$estimates$estimate)))
  expect_true(all(is.na(fit$estimates$std_error)))
  expect_true(all(is.na(fit$covariance$aipw_naive)))
  expect_true(all(is.na(summary(fit)[["95% CI"]])))
  expect_output(print(fit), "cross-fitting used")
  expect_output(print(fit), "Point estimates only")
  expect_true(all(fit$dml$dat$S == 1 & fit$dml$dat$R == 1))

  # Changes outside the responder sample cannot change the fit.
  d[d$S == 0 | d$R == 0, c("X", "W", "Y")] <- NA_real_
  expect_warning(changed <- run(d), "point estimates only")
  expect_equal(fit, changed)
})

test_that("naive XGBoost follows the responder Hajek AIPW formula", {
  expect_warning(fit <- .naive_xgboost(example_data(), c("X", "W"),
    5L, list(nrounds = 3L), 1L), "point estimates only")
  risks <- matrix(0, 5L, 2L)
  for (k in 1:5) {
    d <- fit$dat[fit$dat$fold == k, ]
    for (a in 0:1) {
      mu <- d[[paste0("muhat_", a)]]
      obs <- d$A == a & d$C == 0
      risks[k, a + 1] <- mean(mu) +
        sum((d$Y[obs] - mu[obs]) / d$pihat[obs]) / sum(1 / d$pihat[obs])
    }
  }
  expect_equal(unname(fit$eta_hat), colMeans(risks))
  colnames(risks) <- c("etahat_0", "etahat_1")
  expect_equal(as.matrix(fit$fold_estimates), risks)
})

test_that("naive nuisance fits exclude held-out rows and auxiliary data", {
  calls <- list()
  local_mocked_bindings(
    nonpar_est = function(x, y, method, arguments, wts = NULL) {
      calls[[length(calls) + 1L]] <<- rownames(x)
      NULL
    },
    nonpar_pred = function(mod, newdata, method) rep(.5, nrow(newdata))
  )
  expect_warning(fit <- .naive_xgboost(example_data(), c("X", "W"),
    2L, NULL, 1L), "point estimates only")
  d <- fit$dat
  expect_length(calls, 6L)
  for (k in 1:2) {
    training <- d$fold != k
    j <- 3L * (k - 1L)
    expect_equal(calls[[j + 1L]], rownames(d)[training & d$A == 0 & d$C == 0])
    expect_equal(calls[[j + 2L]], rownames(d)[training & d$A == 1 & d$C == 0])
    expect_equal(calls[[j + 3L]], rownames(d)[training])
  }
})

test_that("naive held-out outcomes cannot affect their own predictions", {
  d <- example_data()
  expect_warning(fit <- .naive_xgboost(d, c("X", "W"), 2L,
    list(nrounds = 3L), 1L), "point estimates only")
  held_out <- fit$dat$fold == 1L
  rows <- rownames(fit$dat)[held_out & fit$dat$C == 0L]
  d[rows, "Y"] <- 1 - d[rows, "Y"]
  expect_warning(changed <- .naive_xgboost(d, c("X", "W"), 2L,
    list(nrounds = 3L), 1L), "point estimates only")
  columns <- c("muhat_0", "muhat_1", "pihat")
  expect_equal(fit$dat[held_out, columns], changed$dat[held_out, columns])
})

test_that("naive cross-fitting rejects folds without enough observed outcomes", {
  d <- example_data()
  expect_error(.naive_xgboost(d, c("X", "W"), 1L, NULL, 1L), "at least 2")
  expect_error(.naive_xgboost(d, c("X", "W"), 2.5, NULL, 1L), "integer")
  expect_error(.naive_xgboost(d, c("X", "W"), 1000L, NULL, 1L),
               "at least K uncensored responders")
})

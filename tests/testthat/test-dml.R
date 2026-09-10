test_that("folds are stratified and reproducible", {
  d <- example_data()
  first <- dml_fit(d, c("X", "W"), c("X", "W"), arguments = list(nrounds = 3L))
  second <- dml_fit(d, c("X", "W"), c("X", "W"), arguments = list(nrounds = 3L))
  expect_true(all(table(first$dat$S, first$dat$A, first$dat$fold) == 50L))
  expect_equal(first, second)
})

test_that("each fold fits the four intended training samples", {
  calls <- list()
  local_mocked_bindings(
    nonpar_est = function(x, y, method, arguments, wts = NULL) {
      calls[[length(calls) + 1L]] <<- as.integer(rownames(x))
      NULL
    },
    nonpar_pred = function(mod, newdata, method) rep(.5, nrow(newdata))
  )
  fit <- dml_fit(example_data(), c("X", "W"), c("X", "W"))
  d <- fit$dat
  expect_length(calls, 20L)
  for (k in 1:5) {
    training <- d$fold != k
    control <- d$Q == 1 & d$A == 0
    treated <- d$Q == 1 & d$A == 1
    j <- 4 * (k - 1)
    expect_equal(calls[[j + 1]], which(training & control))
    expect_equal(calls[[j + 2]], which(training & treated))
    expect_equal(calls[[j + 3]], which(training & (control | d$S == 0)))
    expect_equal(calls[[j + 4]], which(training & (treated | d$S == 0)))
  }
})

test_that("DML uses fold-specific weighted Hajek sums and individual covariance", {
  d <- example_data()
  d$wt[d$S == 0] <- seq(.5, 1.5, length.out = sum(d$S == 0))
  result <- dml_fit(d, c("X", "W"), c("X", "W"), arguments = list(nrounds = 3L))
  fold_risks <- matrix(0, 5, 2)
  covariance <- matrix(0, 2, 2)
  for (k in 1:5) {
    x <- result$dat[result$dat$fold == k, ]
    contribution <- matrix(0, nrow(x), 2L)
    for (a in 0:1) {
      obs <- x$Q == 1 & x$A == a
      aux <- x$S == 0
      mu <- x[[paste0("muhat_", a)]]
      contribution[obs, a + 1] <- (x$Y[obs] - mu[obs]) / x$pihat[obs] / sum(1 / x$pihat[obs])
      contribution[aux, a + 1] <- x$wt[aux] * mu[aux] / sum(x$wt[aux])
    }
    fold_risks[k, ] <- colSums(contribution)
    covariance <- covariance + cov(nrow(x) * contribution) / nrow(d)
  }
  expect_equal(unname(result$eta_hat), colMeans(fold_risks))
  expect_equal(result$eta_hat_cov, covariance / 5)
})

test_that("held-out outcomes do not train their own nuisance models", {
  d <- example_data()
  fit <- dml_fit(d, c("X", "W"), c("X", "W"), arguments = list(nrounds = 3L))
  held_out <- fit$dat$fold == 1L
  observed <- held_out & fit$dat$Q == 1
  d$Y[observed] <- 1 - d$Y[observed]
  changed <- dml_fit(d, c("X", "W"), c("X", "W"), arguments = list(nrounds = 3L))
  columns <- c("muhat_0", "muhat_1", "Q_prob")
  expect_equal(fit$dat[held_out, columns], changed$dat[held_out, columns])
})

test_that("DML does not use unavailable nonresponder covariates", {
  d <- example_data()
  d[d$S == 1 & d$R == 0, c("X", "W")] <- NA_real_
  result <- dml_fit(d, c("X", "W"), c("X", "W"), arguments = list(nrounds = 3L))
  expect_true(all(is.finite(result$eta_hat)))
  expect_true(all(is.finite(result$eta_hat_cov)))
})

test_that("SuperLearner remains an optional DML learner", {
  skip_if_not_installed("SuperLearner")
  result <- dml_fit(example_data(), c("X", "W"), c("X", "W"), K = 2L,
    method = "SuperLearner", arguments = list(SL.library = "SL.glm", cvControl = list(V = 2L)))
  expect_true(all(is.finite(result$eta_hat)))
})

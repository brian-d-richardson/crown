test_that("proposed AIPW keeps the Hajek point estimate and has finite covariance", {
  set.seed(20260915)
  n_trial <- 120L
  n_aux <- 100L
  trial <- data.frame(
    S = 1L, R = 1L, C = rbinom(n_trial, 1, 0.2),
    A = rbinom(n_trial, 1, 0.5), W = rnorm(n_trial), wt = 1
  )
  trial$Y <- rbinom(
    n_trial, 1, plogis(-0.4 + 0.7 * trial$A + 0.3 * trial$W)
  )
  trial$Y[trial$C == 1L] <- 0L
  auxiliary <- data.frame(
    S = 0L, R = 0L, C = 0L, A = rbinom(n_aux, 1, 0.5),
    W = rnorm(n_aux, 0.3), Y = 0L, wt = runif(n_aux, 0.5, 1.5)
  )
  dat <- rbind(trial, auxiliary)

  result <- fit_aipw(dat, Y ~ A * W, Q ~ W, "proposed")

  dat$Q <- dat$S * dat$R * (1L - dat$C)
  observed <- dat$Q == 1L
  auxiliary_rows <- dat$S == 0L
  outcome_fit <- glm(Y ~ A * W, binomial(), dat[observed, ])
  dat0 <- dat1 <- dat
  dat0$A <- 0L
  dat1$A <- 1L
  mu0 <- predict(outcome_fit, dat0, type = "response")
  mu1 <- predict(outcome_fit, dat1, type = "response")

  manual <- vapply(0:1, function(a) {
    observed_a <- observed & dat$A == a
    selection_fit <- glm(
      Q ~ W, binomial(), dat[observed_a | auxiliary_rows, ],
      weights = dat$wt[observed_a | auxiliary_rows]
    )
    probability <- predict(selection_fit, dat[observed_a, ], type = "response")
    odds <- probability / (1 - probability)
    mu <- if (a == 0L) mu0 else mu1
    weighted.mean(mu[auxiliary_rows], dat$wt[auxiliary_rows]) +
      weighted.mean(dat$Y[observed_a] - mu[observed_a], 1 / odds)
  }, numeric(1L))

  expect_equal(unname(unlist(result[c("etahat_0", "etahat_1")])), manual)
  covariance <- matrix(
    c(result$cov_00, result$cov_01, result$cov_01, result$cov_11), 2L
  )
  expect_true(all(is.finite(covariance)))
  expect_true(all(diag(covariance) > 0))
  expect_true(all(eigen(covariance, symmetric = TRUE)$values > 0))
})

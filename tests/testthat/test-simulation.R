test_that("Simulation 2 runs all supported estimators and agrees with DML", {
  path <- test_path("..", "..", "simulation", "sim_functions", "sim2.R")
  skip_if_not(file.exists(path), "Simulation scripts are source-only")
  env <- new.env(parent = asNamespace("crown"))
  sys.source(path, env)
  result <- env$run_sim2(
    data.frame(n_trial = 500L, n_auxiliary = 500L),
    "test", tempdir(), mc_reps = 2L, arguments = list(nrounds = 3L)
  )
  expect_equal(nrow(result$results), 14L)
  expect_true(file.exists(result$result_file))
  expect_true(all(is.finite(result$results$rdhat)))
  expect_equal(
    unique(result$results[c("Version", "Model", "Estimator")]),
    data.frame(
      Version = c(rep("Proposed", 3L), rep("Naive", 3L), "Proposed"),
      Model = c(rep("Parametric", 6L), "Nonparametric"),
      Estimator = c(
        "G-Formula", "IPW", "AIPW", "G-Formula", "IPW", "AIPW", "AIPW"
      )
    )
  )
  d <- env$generate_sim2_data(20L, 500L, 500L, .5, .3, 91000001L)$data
  fit <- dml_fit(d, c("X1", "W1", "W2", "W3"), c("X1", "W1", "W2", "W3"),
    arguments = list(nrounds = 3L), random_seed = 91000001L)
  nonparametric <- subset(
    result$results,
    seed == 91000001L & Model == "Nonparametric"
  )
  expect_equal(nonparametric$etahat_0, unname(fit$eta_hat[1]))
  expect_equal(nonparametric$etahat_1, unname(fit$eta_hat[2]))
})

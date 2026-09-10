test_that("Simulation 2 runs sequentially and agrees with a direct fit", {
  path <- test_path("..", "..", "simulation", "sim_functions", "sim2.R")
  skip_if_not(file.exists(path), "Simulation scripts are source-only")
  env <- new.env(parent = asNamespace("crown"))
  sys.source(path, env)
  result <- env$run_sim2(
    data.frame(n_trial = 500L, n_auxiliary = 500L),
    "test", tempdir(), mc_reps = 2L, arguments = list(nrounds = 3L)
  )
  expect_equal(nrow(result$results), 4L)
  expect_true(file.exists(result$result_file))
  expect_true(all(is.finite(result$results$rdhat)))
  d <- env$generate_sim2_data(20L, 500L, 500L, .5, .3, 91000001L)$data
  fit <- dml_fit(d, c("X1", "W1", "W2", "W3"), c("X1", "W1", "W2", "W3"),
    arguments = list(nrounds = 3L), random_seed = 91000001L)
  expect_equal(result$results$etahat_0[2], unname(fit$eta_hat[1]))
  expect_equal(result$results$etahat_1[2], unname(fit$eta_hat[2]))
})

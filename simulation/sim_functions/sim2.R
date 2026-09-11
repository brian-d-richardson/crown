# Simulation 2: all supported estimators ------------------------------------

# Generate one trial and one auxiliary sample from the nonlinear DGP.
generate_sim2_data <- function(
    n_clusters, n_trial, n_auxiliary, p_response, p_censoring, seed) {

  stopifnot(n_clusters == 20L, p_response == 0.5, p_censoring == 0.3)
  n_trial <- n_clusters * round(n_trial / n_clusters)
  n_auxiliary <- n_clusters * round(n_auxiliary / n_clusters)
  set.seed(seed)

  cluster_risk <- c(
    0.92, 0.78, 0.07, -1.99, 0.62, -0.06, -0.16, -1.47, -0.48, 0.42,
    1.36, -0.10, 0.39, -0.05, -1.38, -0.41, -0.39, -0.06, 1.10, 0.76
  )
  treatment_by_cluster <- sample(rep(c(0, 1), n_clusters / 2))
  response_intercept <- -1.417151
  censoring_intercept <- -0.8532847

  cluster <- rep(seq_len(n_clusters), length.out = n_trial)
  A <- treatment_by_cluster[cluster]
  X1 <- cluster_risk[cluster]
  W1 <- rbinom(n_trial, 1, 0.5)
  W2 <- rnorm(n_trial)
  W3 <- rnorm(n_trial)
  R <- rbinom(
    n_trial, 1,
    plogis(response_intercept + 2 * (W2^2 < 1))
  )
  C <- rbinom(
    n_trial, 1,
    plogis(censoring_intercept - 0.25 * A + 0.25 * W1)
  )
  C[R == 0L] <- 0L
  mu0 <- ifelse(W2^2 < 1, 0.9, plogis(0.25 * sin(pi * W3 / 4)))
  mu1 <- ifelse(W2^2 < 1, 0.1, plogis(0.25 * sin(pi * W3 / 4)))
  Y0 <- rbinom(n_trial, 1, mu0)
  Y1 <- rbinom(n_trial, 1, mu1)
  Y <- (1 - A) * Y0 + A * Y1
  Y[R == 0L | C == 1L] <- 0L

  trial <- data.frame(
    id = seq_len(n_trial), cluster, X1, W1, W2, W3, A, R, C, Y, wt = 1,
    S = 1L
  )

  cluster <- rep(seq_len(n_clusters), length.out = n_auxiliary)
  W1 <- rbinom(n_auxiliary, 1, 0.75)
  auxiliary <- data.frame(
    id = n_trial + seq_len(n_auxiliary),
    cluster,
    X1 = cluster_risk[cluster],
    W1,
    W2 = rnorm(n_auxiliary),
    W3 = rnorm(n_auxiliary),
    A = treatment_by_cluster[cluster],
    R = 0L,
    C = 0L,
    Y = 0L,
    wt = ifelse(W1 == 0L, 0.75, 0.25),
    S = 0L
  )
  auxiliary$wt <- auxiliary$wt / mean(auxiliary$wt)

  list(
    data = rbind(trial, auxiliary),
    truth = c(eta0 = mean(mu0), eta1 = mean(mu1)),
    n_trial = n_trial,
    n_auxiliary = n_auxiliary
  )
}

# Run the Monte Carlo replicates.
run_sim2 <- function(
    sample_sizes, run_id, out_dir, mc_reps = 20L, base_seed = 91000000L,
    K = 5L, arguments = NULL) {

  # Monte Carlo settings.
  grid <- sample_sizes[rep(seq_len(nrow(sample_sizes)), mc_reps), ]
  grid$replicate <- rep(seq_len(mc_reps), each = nrow(sample_sizes))
  grid$seed <- base_seed + seq_len(nrow(grid))
  covariates <- c("X1", "W1", "W2", "W3")
  results <- data.frame()
  started_at <- proc.time()[["elapsed"]]

  for (i in seq_len(nrow(grid))) {

    # Generate one trial and one auxiliary sample.
    generated <- generate_sim2_data(
      20L, grid$n_trial[i], grid$n_auxiliary[i], 0.5, 0.3, grid$seed[i]
    )

    # Fit the six parametric combinations.
    outcome_formula <- Y ~ A * (X1 + W1 + W2 + W3)
    propensity_formula <- Q ~ X1 + W1 + W2 + W3
    censoring_formula <- C ~ X1 + W1 + W2 + W3
    parametric <- rbind(
      fit_gformula(generated$data, outcome_formula, "proposed"),
      fit_ipw(generated$data, propensity_formula, "proposed"),
      fit_aipw(
        generated$data, outcome_formula, propensity_formula, "proposed"
      ),
      fit_gformula(generated$data, outcome_formula, "naive"),
      fit_ipw(generated$data, censoring_formula, "naive"),
      fit_aipw(
        generated$data, outcome_formula, censoring_formula, "naive"
      )
    )

    # Fit the supported Proposed nonparametric AIPW estimator.
    dml <- dml_fit(
      generated$data, covariates, covariates, K,
      arguments, random_seed = grid$seed[i]
    )

    # Store risks, contrasts, and estimated variances.
    result <- rbind(
      parametric,
      .eta_result(dml$eta_hat[1], dml$eta_hat[2], dml$eta_hat_cov)
    )
    result$Estimator <- c(
      "G-Formula", "IPW", "AIPW",
      "G-Formula", "IPW", "AIPW", "AIPW"
    )
    result$Version <- c(rep("Proposed", 3L), rep("Naive", 3L), "Proposed")
    result$Model <- c(rep("Parametric", 6L), "Nonparametric")
    result$rdhat <- result$etahat_1 - result$etahat_0
    result$rrhat <- result$etahat_1 / result$etahat_0
    result$var_rd <- result$cov_00 + result$cov_11 - 2 * result$cov_01
    result$var_rr <- result$cov_11 / result$etahat_0^2 +
      result$cov_00 * result$etahat_1^2 / result$etahat_0^4 -
      2 * result$cov_01 * result$etahat_1 / result$etahat_0^3

    result$eta_0 <- generated$truth[1]
    result$eta_1 <- generated$truth[2]
    result$rd <- generated$truth[2] - generated$truth[1]
    result$rr <- generated$truth[2] / generated$truth[1]
    result$seed <- grid$seed[i]
    result$n_trial <- grid$n_trial[i]
    result$n_aux <- grid$n_auxiliary[i]
    result$K <- K
    result$run_id <- run_id
    result$replicate <- grid$replicate[i]
    results <- rbind(results, result)
    cat("Completed", i, "of", nrow(grid), "replicates\n")
  }

  # Save results and total running time.
  elapsed <- proc.time()[["elapsed"]] - started_at
  rownames(results) <- NULL
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(out_dir, paste0("monte_carlo_", run_id, "_results.csv"))
  write.csv(results, result_file, row.names = FALSE)
  timing <- data.frame(mc_reps, combinations = nrow(sample_sizes), K,
                       elapsed_seconds = elapsed)
  write.csv(timing, file.path(out_dir, paste0(run_id, "_timing.csv")), row.names = FALSE)
  cat("Finished in", round(elapsed, 1), "seconds\n")
  invisible(list(results = results, result_file = result_file, elapsed = elapsed))
}

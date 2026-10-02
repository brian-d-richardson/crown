# Simulation 1: model misspecification --------------------------------------

# Generate one trial and one auxiliary sample.
generate_sim1_data <- function(
    n_clusters, n_trial, n_auxiliary, p_response, p_censoring, seed) {

  response_intercept <- uniroot(
    function(value) {
      mean(plogis(c(value, value - 2, value, value))) - p_response
    },
    c(-20, 20)
  )$root
  combinations <- expand.grid(A = c(0, 1), W1 = c(0, 1))
  response_probability <- plogis(
    response_intercept + combinations$A - combinations$W1
  )
  censoring_intercept <- uniroot(
    function(value) {
      censoring_probability <- plogis(
        value - 0.25 * combinations$A + 0.25 * combinations$W1
      )
      sum(censoring_probability * response_probability * 0.25) /
        sum(response_probability * 0.25) - p_censoring
    },
    c(-20, 20)
  )$root

  n_trial <- n_clusters * round(n_trial / n_clusters)
  n_auxiliary <- n_clusters * round(n_auxiliary / n_clusters)
  set.seed(seed)
  cluster_risk <- c(
    0.92, 0.78, 0.07, -1.99, 0.62, -0.06, -0.16, -1.47, -0.48, 0.42,
    1.36, -0.10, 0.39, -0.05, -1.38, -0.41, -0.39, -0.06, 1.10, 0.76
  )
  treatment_by_cluster <- sample(rep(c(0L, 1L), n_clusters / 2L))

  cluster <- rep(seq_len(n_clusters), length.out = n_trial)
  A <- treatment_by_cluster[cluster]
  X1 <- cluster_risk[cluster]
  W1 <- rbinom(n_trial, 1, 0.5)
  W2 <- rnorm(n_trial)
  R <- rbinom(
    n_trial, 1, plogis(response_intercept - 2 * A * W1)
  )
  C <- rbinom(
    n_trial, 1, plogis(censoring_intercept - 0.25 * A + 0.25 * W1)
  )
  C[R == 0L] <- 0L
  probability0 <- plogis(-1 + 2 * W1 + 0.5 * W2 + 0.25 * X1)
  probability1 <- plogis(-W1 - 0.5 * W2)
  Y0 <- rbinom(n_trial, 1, probability0)
  Y1 <- rbinom(n_trial, 1, probability1)
  Y <- (1 - A) * Y0 + A * Y1
  Y[R == 0L | C == 1L] <- 0L
  trial <- data.frame(
    id = seq_len(n_trial), cluster, X1, W1, W2, A, R, C, Y, wt = 1, S = 1L
  )

  cluster <- rep(seq_len(n_clusters), length.out = n_auxiliary)
  A <- treatment_by_cluster[cluster]
  W1 <- rbinom(n_auxiliary, 1, 0.75)
  auxiliary <- data.frame(
    id = n_trial + seq_len(n_auxiliary),
    cluster,
    X1 = cluster_risk[cluster],
    W1,
    W2 = rnorm(n_auxiliary),
    A,
    R = 0L,
    C = 0L,
    Y = 0L,
    wt = ifelse(W1 == 0L, 0.75, 0.25),
    S = 0L
  )

  list(
    data = rbind(trial, auxiliary),
    truth = c(eta0 = mean(Y0), eta1 = mean(Y1)),
    n_trial = n_trial,
    n_auxiliary = n_auxiliary
  )
}

# Fit the three proposed parametric estimators for one model scenario.
fit_sim1_estimators <- function(data, mu_correct, pi_correct) {
  outcome_formula <- if (mu_correct) {
    Y ~ A * (X1 + W1 + W2)
  } else {
    Y ~ A * (X1 + W2)
  }
  propensity_formula <- if (pi_correct) {
    Q ~ X1 + W1 + W2
  } else {
    Q ~ X1 + W2
  }

  result <- rbind(
    fit_gformula(data, outcome_formula, "proposed"),
    fit_ipw(data, propensity_formula, "proposed"),
    fit_aipw(data, outcome_formula, propensity_formula, "proposed")
  )
  result$Estimator <- c("G-Formula", "IPW", "AIPW")
  result$Version <- "Proposed"
  result
}

# Add truths, contrasts, and variance estimates to fitted risks.
.complete_sim1_result <- function(
    result, generated, seed, n_clusters, p_response, p_censoring,
    mu_correct, pi_correct, K) {
  result$rdhat <- result$etahat_1 - result$etahat_0
  result$rrhat <- result$etahat_1 / result$etahat_0
  result$var_rd <- result$cov_00 + result$cov_11 - 2 * result$cov_01
  result$var_rr <- result$cov_11 / result$etahat_0^2 +
    result$cov_00 * result$etahat_1^2 / result$etahat_0^4 -
    2 * result$cov_01 * result$etahat_1 / result$etahat_0^3
  result$eta_0 <- generated$truth[["eta0"]]
  result$eta_1 <- generated$truth[["eta1"]]
  result$rd <- result$eta_1 - result$eta_0
  result$rr <- result$eta_1 / result$eta_0
  result$seed <- seed
  result$m <- n_clusters
  result$n_trial <- generated$n_trial
  result$n_aux <- generated$n_auxiliary
  result$p_resp <- p_response
  result$p_cens <- p_censoring
  result$mu_correct <- mu_correct
  result$pi_correct <- pi_correct
  result$K <- K
  result
}

# Run the Monte Carlo study.
run_sim1 <- function(
    sample_sizes, run_id, out_dir,
    n_clusters = 20L, p_response = 0.5, p_censoring = 0.3,
    mc_reps = 30L, base_seed = 11000000L, replicates = seq_len(mc_reps),
    K = 5L, arguments = NULL) {

  # Monte Carlo settings.

  grid <- sample_sizes[rep(seq_len(nrow(sample_sizes)), length(replicates)), ]
  grid$replicate <- rep(replicates, each = nrow(sample_sizes))
  grid$seed <- base_seed + grid$replicate
  rownames(grid) <- NULL
  scenarios <- expand.grid(
    mu_correct = c(TRUE, FALSE),
    pi_correct = c(TRUE, FALSE)
  )

  cat(
    "Running", length(replicates), "Simulation 1 replicates for each of",
    nrow(sample_sizes), "sample-size combinations\n"
  )

  # Run the simulations.

  results <- data.frame()
  for (i in seq_len(nrow(grid))) {
    generated <- generate_sim1_data(
      n_clusters, grid$n_trial[i], grid$n_auxiliary[i],
      p_response, p_censoring, grid$seed[i]
    )

    for (j in seq_len(nrow(scenarios))) {
      estimates <- fit_sim1_estimators(
        generated$data, scenarios$mu_correct[j], scenarios$pi_correct[j]
      )
      results <- rbind(results, .complete_sim1_result(
        estimates, generated, grid$seed[i], n_clusters, p_response, p_censoring,
        scenarios$mu_correct[j], scenarios$pi_correct[j], K
      ))
    }

    correct_outcome <- Y ~ A * (X1 + W1 + W2)
    correct_censoring <- C ~ X1 + W1 + W2
    naive <- rbind(
      fit_gformula(generated$data, correct_outcome, "naive"),
      fit_ipw(generated$data, correct_censoring, "naive"),
      fit_aipw(generated$data, correct_outcome, correct_censoring, "naive")
    )
    naive$Estimator <- c("G-Formula", "IPW", "AIPW")
    naive$Version <- "Naive"
    results <- rbind(results, .complete_sim1_result(
      naive, generated, grid$seed[i], n_clusters, p_response, p_censoring,
      TRUE, TRUE, K
    ))

    dml <- dml_fit(
      generated$data, c("X1", "W1", "W2"), c("X1", "W1", "W2"),
      K, arguments, random_seed = grid$seed[i]
    )
    dml <- .eta_result(dml$eta_hat[1], dml$eta_hat[2], dml$eta_hat_cov)
    dml$Estimator <- "DML"
    dml$Version <- "Proposed"
    results <- rbind(results, .complete_sim1_result(
      dml, generated, grid$seed[i], n_clusters, p_response, p_censoring,
      TRUE, TRUE, K
    ))
  }

  # Save the results.

  rownames(results) <- NULL
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(out_dir, paste0(run_id, ".csv"))
  write.csv(results, result_file, row.names = FALSE)
  invisible(list(results = results, result_file = result_file))
}

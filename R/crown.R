#' Cluster-Randomized Trial Analysis with Outcomes Weighted for Nonresponse
#'
#' G-formula, IPW, and AIPW estimators combine trial outcomes with auxiliary
#' covariates using logistic regression or cross-fitted machine learning.
#' @importFrom SuperLearner SuperLearner
#' @importFrom numDeriv jacobian
#' @importFrom stats as.formula binomial coef cov glm model.frame model.matrix model.response plogis predict weighted.mean
#' @importFrom utils tail
#' @importFrom xgboost xgboost
#' @keywords internal
"_PACKAGE"

#' Analyze a trial and an auxiliary sample
#'
#' Runs one selected estimator using logistic regression or cross-fitted XGBoost.
#'
#' @param trial_data,auxiliary_data Data frames for the trial and auxiliary sample.
#' @param outcome Binary outcome column; required for uncensored responders.
#' @param treatment Cluster-level treatment column in the trial. Auxiliary
#'   treatment is matched using the cluster identifier.
#' @param response Trial response indicator column (1 = responder).
#' @param censoring Censoring indicator column (1 = censored).
#' @param cluster Cluster identifier column shared by both samples. Treatment
#'   must be constant within each cluster.
#' @param covariates Shared baseline covariate column names, including the
#'   cluster-level covariates needed for conditional independence.
#' @param version `"proposed"` (default) or `"naive"`. Proposed estimators use
#'   auxiliary covariates; naive estimators target trial responders.
#' @param estimator `"aipw"` (default), `"gformula"`, or `"ipw"`.
#'   Only the selected estimator and its required nuisance models are fitted.
#' @param model `"xgboost"` (default) uses cross-fitting; `"logistic"` uses
#'   logistic regressions. Only Proposed AIPW has XGBoost standard errors.
#' @param auxiliary_weight Optional auxiliary sampling-weight column. Omit for
#'   equal weights. Supplied weights are normalized to have mean one.
#' @param treatment_values Control and treatment values, in that order.
#' @param outcome_formula,propensity_formula,censoring_formula Optional logistic
#'   formulas using standardized names `Y`, `A`, `Q`, `C`, and the covariates.
#'   Defaults are `Y ~ A * (covariates)`, `Q ~ covariates`, and `C ~ covariates`.
#' @param K Number of outer cross-fitting folds (default 5); used only by
#'   XGBoost, not tuning CV.
#' @param arguments Optional XGBoost settings, e.g., `list(nrounds = 100L)`;
#'   used only by XGBoost.
#' @param random_seed Seed set before creating folds and fitting XGBoost.
#'
#' @details Logistic estimators use sandwich covariance. Proposed AIPW with
#'   XGBoost uses cross-fitting and provides standard errors. Other XGBoost
#'   combinations return point estimates only.
#'
#' @return A `crown_fit` list containing:
#' \itemize{
#'   \item `estimates`: risks under control and treatment, RD, RR, standard
#'     errors, and 95 percent Wald confidence intervals (standard errors and
#'     intervals are `NA` for XGBoost except Proposed AIPW);
#'   \item `variance`: variances on the same parameter scales;
#'   \item `covariance`: one two-by-two risk covariance matrix per estimator;
#'   \item `dml`: fold-specific estimates and out-of-fold predictions, included
#'     only for `model = "xgboost"`.
#' }
#' @export
crown <- function(
    trial_data,
    auxiliary_data,
    outcome,
    treatment,
    response,
    censoring,
    cluster,
    covariates,
    version = c("proposed", "naive"),
    estimator = c("aipw", "gformula", "ipw"),
    model = c("xgboost", "logistic"),
    auxiliary_weight = NULL,
    treatment_values = c(0, 1),
    outcome_formula = NULL,
    propensity_formula = NULL,
    censoring_formula = NULL,
    K = 5L,
    arguments = NULL,
    random_seed = 1L) {

  version <- match.arg(version)
  model <- match.arg(model)
  estimator <- match.arg(estimator)

  # Prepare the combined trial and auxiliary data.
  data <- .prepare_crown_data(
    trial_data, auxiliary_data, outcome, treatment, response, censoring,
    cluster, covariates, auxiliary_weight, treatment_values
  )

  dml <- NULL
  if (model == "logistic") {

    # Specify the logistic outcome and weighting models.
    covariate_terms <- paste(sprintf("`%s`", covariates), collapse = " + ")
    if (is.null(outcome_formula)) {
      outcome_formula <- as.formula(paste0("Y ~ A * (", covariate_terms, ")"))
    }
    if (is.null(propensity_formula)) {
      propensity_formula <- as.formula(paste("Q ~", covariate_terms))
    }
    if (is.null(censoring_formula)) {
      censoring_formula <- as.formula(paste("C ~", covariate_terms))
    }
    weight_formula <- if (version == "naive") censoring_formula else propensity_formula

    # Fit only the selected estimator.
    if (estimator == "gformula") {
      result <- fit_gformula(data, outcome_formula, version)
    } else if (estimator == "ipw") {
      result <- fit_ipw(data, weight_formula, version)
    } else {
      result <- fit_aipw(data, outcome_formula, weight_formula, version)
    }

  } else {

    # Fit the selected cross-fitted estimator.
    if (version == "proposed") {
      dml <- dml_fit(
        data, covariates, covariates, K, "xgboost", arguments, random_seed,
        estimator = estimator
      )
    } else {
      dml <- .naive_xgboost(data, covariates, K, arguments, random_seed, estimator)
    }
    result <- .eta_result(
      dml$eta_hat[[1]], dml$eta_hat[[2]], dml$eta_hat_cov
    )
  }

  # Return the estimates, standard errors, and confidence intervals.
  results <- list(result)
  names(results) <- paste(estimator, version, sep = "_")
  output <- .format_crown_results(results)
  output$dml <- dml
  class(output) <- c("crown_fit", "list")
  output
}

#' Standardize trial and auxiliary columns
#' @inheritParams crown
#' @return Combined data with standardized indicator names and normalized weights.
#' @keywords internal
#' @noRd
.prepare_crown_data <- function(
    trial_data, auxiliary_data, outcome, treatment, response, censoring,
    cluster, covariates, auxiliary_weight, treatment_values) {
  A <- match(trial_data[[treatment]], treatment_values) - 1L
  assignment <- unique(data.frame(cluster = trial_data[[cluster]], A = A))
  if (anyNA(A) || anyDuplicated(assignment$cluster)) {
    stop("Each trial cluster must have one treatment value from treatment_values.")
  }
  aux_A <- A[match(auxiliary_data[[cluster]], trial_data[[cluster]])]
  if (anyNA(aux_A)) stop("Auxiliary cluster identifiers must match trial clusters.")

  R <- trial_data[[response]]
  C <- trial_data[[censoring]]
  Y <- trial_data[[outcome]]
  C[R == 0L] <- 0L
  Y[R == 0L | C == 1L] <- 0L

  auxiliary_weights <- if (is.null(auxiliary_weight)) {
    rep(1, nrow(auxiliary_data))
  } else {
    auxiliary_data[[auxiliary_weight]]
  }
  auxiliary_weights <- auxiliary_weights / mean(auxiliary_weights)

  trial <- cbind(
    trial_data[covariates],
    data.frame(.cluster = trial_data[[cluster]], A = A, S = 1L,
               R = R, C = C, Y = Y, wt = 1)
  )
  auxiliary <- cbind(
    auxiliary_data[covariates],
    data.frame(.cluster = auxiliary_data[[cluster]], A = aux_A, S = 0L,
               R = 0L, C = 0L, Y = 0L, wt = auxiliary_weights)
  )
  rbind(trial, auxiliary)
}

#' Format causal estimates and delta-method variances
#' @param results Named list of rows from the individual estimators.
#' @return Estimate and variance tables, and named risk covariance matrices.
#' @keywords internal
#' @noRd
.format_crown_results <- function(results) {
  output <- list(estimates = data.frame(), variance = data.frame(), covariance = list())
  for (name in names(results)) {

    # Arm-specific risks and covariance.
    row <- results[[name]]
    eta <- c(row$etahat_0, row$etahat_1)
    covariance <- matrix(c(row$cov_00, row$cov_01, row$cov_01, row$cov_11), 2)

    # Risk difference, risk ratio, and delta-method variances.
    estimate <- c(
      mean_control = eta[1],
      mean_treated = eta[2],
      risk_difference = eta[2] - eta[1],
      risk_ratio = eta[2] / eta[1]
    )
    gradient <- rbind(
      c(1, 0), c(0, 1), c(-1, 1),
      c(-eta[2] / eta[1]^2, 1 / eta[1])
    )
    variance <- diag(gradient %*% covariance %*% t(gradient))
    standard_error <- sqrt(variance)

    # Estimator and version labels.
    label <- strsplit(name, "_", fixed = TRUE)[[1]]
    estimator <- c(gformula = "G-Formula", ipw = "IPW", aipw = "AIPW", dml = "AIPW")[[label[1]]]
    version <- c(naive = "Naive", proposed = "Proposed")[[label[2]]]
    estimates <- data.frame(
        estimator, version, parameter = names(estimate),
        estimate = unname(estimate), std_error = standard_error,
        conf_low = unname(estimate) - 1.96 * standard_error,
        conf_high = unname(estimate) + 1.96 * standard_error
    )
    variances <- data.frame(estimator, version, parameter = names(estimate), variance)
    output$estimates <- rbind(output$estimates, estimates)
    output$variance <- rbind(output$variance, variances)
    output$covariance[[name]] <- covariance
  }
  output
}

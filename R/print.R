#' Print a Crown fit
#'
#' Displays the nuisance-model type and four causal estimates for each
#' estimator and version. XGBoost fits also report the number of cross-fitting
#' folds.
#'
#' @param x A `crown_fit` object returned by [crown()].
#' @param digits Number of significant digits used to format the printed table.
#' @param ... Additional arguments; currently unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.crown_fit <- function(x, digits = 3, ...) {
  estimates <- x$estimates
  if (x$model == "xgboost") {
    estimates$estimator[
      estimates$estimator == "AIPW" & estimates$version == "Proposed"
    ] <- "DML"
  }
  estimates$model <- if (x$model == "xgboost") "XGBoost" else "Logistic"
  output <- unique(estimates[c("estimator", "version", "model")])
  parameters <- c(
    mean_control = "eta(0)",
    mean_treated = "eta(1)",
    risk_difference = "RD",
    risk_ratio = "RR"
  )
  output_key <- paste(output$estimator, output$version)

  for (parameter in names(parameters)) {
    rows <- estimates$parameter == parameter
    estimate_key <- paste(estimates$estimator[rows], estimates$version[rows])
    output[[parameters[[parameter]]]] <- estimates$estimate[rows][
      match(output_key, estimate_key)
    ]
  }

  cat("Crown causal estimates\n\n")
  print(output, row.names = FALSE, digits = digits)
  if (!is.null(x$dml)) {
    cat("\n", nrow(x$dml$fold_estimates),
        " folds of cross-fitting used.\n", sep = "")
  }
  invisible(x)
}


#' Summarize a Crown fit
#'
#' Returns all estimator-specific causal estimates with standard errors and
#' 95 percent confidence intervals.
#'
#' @param object A `crown_fit` object returned by [crown()].
#' @param digits Number of decimal places for estimates, standard errors, and
#'   confidence intervals.
#' @param ... Additional arguments; currently unused.
#'
#' @return A data frame containing the estimator, version, nuisance-model type,
#'   parameter, estimate, standard error, and 95 percent confidence interval.
#'
#' @export
summary.crown_fit <- function(object, digits = 3, ...) {
  output <- object$estimates
  if (object$model == "xgboost") {
    output$estimator[
      output$estimator == "AIPW" & output$version == "Proposed"
    ] <- "DML"
  }
  output$model <- if (object$model == "xgboost") "XGBoost" else "Logistic"
  output[["95% CI"]] <- sprintf(
    "[%.*f, %.*f]",
    digits, output$conf_low,
    digits, output$conf_high
  )
  output[["95% CI"]][is.na(output$conf_low) | is.na(output$conf_high)] <- NA_character_
  output$estimate <- round(output$estimate, digits)
  output$std_error <- round(output$std_error, digits)
  output <- output[c(
    "estimator", "version", "model", "parameter", "estimate", "std_error", "95% CI"
  )]
  rownames(output) <- NULL
  output
}

# Simulation 2 analysis and figures ---------------------------------------

# Summarize and plot the simulation results.
plot_sim2 <- function(results, run_id, figure_dir, dpi = 600) {
  summary <- select(results,
    seed, Version, Estimator, n_trial, n_aux,
    eta_0, eta_1, rd, rr, etahat_0, etahat_1, rdhat, rrhat,
    cov_00, cov_11, var_rd, var_rr
  )
  summary <- rename(summary,
    truth_eta0 = eta_0, truth_eta1 = eta_1,
    truth_rd = rd, truth_rr = rr,
    est_eta0 = etahat_0, est_eta1 = etahat_1,
    est_rd = rdhat, est_rr = rrhat,
    Var_eta0 = cov_00, Var_eta1 = cov_11,
    Var_rd = var_rd, Var_rr = var_rr
  )
  summary <- pivot_longer(summary,
    cols = matches("^(truth|est|Var)_"),
    names_to = c("type", "parameter"),
    names_sep = "_",
    values_to = "value"
  )
  summary <- pivot_wider(summary, names_from = type, values_from = value)
  summary <- mutate(summary,
    lower = est - 1.96 * sqrt(Var),
    upper = est + 1.96 * sqrt(Var)
  )
  summary <- group_by(summary, Version, Estimator, n_trial, n_aux, parameter)
  summary <- summarise(summary,
    empirical_variance = var(est),
    estimated_variance = mean(Var),
    bias = mean(est - truth),
    mse = mean((est - truth)^2),
    coverage = mean(truth >= lower & truth <= upper),
    .groups = "drop"
  )
  summary <- mutate(summary,
    Estimator = factor(Estimator, levels = c("AIPW", "DML")),
    parameter = factor(
      parameter,
      levels = c("eta0", "eta1", "rd", "rr"),
      labels = c("eta(0)", "eta(1)", "RD", "RR")
    )
  )

  sample_size <- function(auxiliary, trial, separator = ".") {
    interaction(
      factor(auxiliary, levels = sort(unique(auxiliary))),
      factor(trial, levels = sort(unique(trial))),
      sep = separator,
      drop = TRUE
    )
  }
  size_labels <- function(labels, separator) {
    values <- strsplit(as.character(labels), separator, fixed = TRUE)
    as.expression(lapply(values, function(value) {
      bquote(n^aux == .(value[[1]]) * "," ~~ n^trial == .(value[[2]]))
    }))
  }

  results <- mutate(results,
    Version = factor(Version, levels = c("Naive", "Proposed")),
    Estimator = factor(Estimator, levels = c("AIPW", "DML")),
    Sample_Size = sample_size(n_aux, n_trial)
  )
  colors <- c("#FF6800", "#803E75", "#C10020", "#FFB300")

  estimate_plot <- ggplot(
    results,
    aes(x = Sample_Size, y = rdhat, color = Estimator, fill = Estimator)
  ) +
    geom_boxplot(alpha = 0.5) +
    geom_hline(yintercept = mean(results$rd), linetype = "dashed") +
    facet_nested(. ~ Estimator) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(
      x = "Auxiliary and Trial Sample Sizes",
      y = expression(hat(RD))
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    guides(x = guide_axis_nested())

  variance <- filter(summary, empirical_variance > 0, estimated_variance > 0)
  variance <- mutate(variance, Sample_Size = sample_size(n_aux, n_trial, "_"))
  variance_plot <- ggplot(
    variance,
    aes(
      x = empirical_variance,
      y = estimated_variance,
      color = Sample_Size,
      shape = parameter
    )
  ) +
    geom_abline(linetype = "dashed") +
    geom_point(size = 3) +
    facet_nested(. ~ Estimator) +
    scale_x_continuous(transform = "log10", breaks = c(0.001, 0.01)) +
    scale_y_continuous(transform = "log10", breaks = c(0.001, 0.01)) +
    scale_color_manual(values = colors, labels = function(x) size_labels(x, "_")) +
    scale_shape_discrete(labels = function(x) parse(text = x)) +
    guides(
      color = guide_legend(nrow = 2, order = 1),
      shape = guide_legend(order = 2)
    ) +
    labs(
      x = "Empirical Variance",
      y = "Average Estimated Variance",
      color = "Sample Size",
      shape = "Parameter"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = unit(-5, "pt")
    )

  coverage <- mutate(summary, Sample_Size = sample_size(n_aux, n_trial))
  coverage_plot <- ggplot(
    coverage,
    aes(x = Sample_Size, y = coverage, color = parameter, shape = parameter)
  ) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    facet_nested(. ~ Estimator) +
    scale_color_manual(values = colors, labels = function(x) parse(text = x)) +
    scale_shape_discrete(labels = function(x) parse(text = x)) +
    labs(
      x = "Auxiliary and Trial Sample Sizes",
      y = "Empirical CI Coverage",
      color = "Estimand",
      shape = "Estimand"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    guides(x = guide_axis_nested())

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(summary, file.path(figure_dir, paste0(run_id, "_summary.csv")), row.names = FALSE)
  files <- file.path(
    figure_dir,
    paste0("monte_carlo_", run_id, c("_estimates.png", "_variance.png", "_confidence.png"))
  )
  ggsave(files[[1]], estimate_plot, width = 5, height = 3, dpi = dpi)
  ggsave(files[[2]], variance_plot, width = 5.8, height = 4, dpi = dpi)
  ggsave(files[[3]], coverage_plot, width = 5, height = 3.5, dpi = dpi)
  invisible(files)
}

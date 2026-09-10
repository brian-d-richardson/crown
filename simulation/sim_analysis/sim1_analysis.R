# Simulation 1 estimates, variances, and confidence coverage -------------

# Summarize and plot the simulation results.
plot_sim1 <- function(results, run_id, figure_dir, dpi = 600) {
  summary <- mutate(results,
    Mu = if_else(
      Version == "Naive", "-",
      if_else(mu_correct, "Correct Mu", "Incorrect Mu")
    ),
    Pi = if_else(
      Version == "Naive", "-",
      if_else(pi_correct, "Correct Pi", "Incorrect Pi")
    )
  )
  summary <- filter(summary, Version == "Proposed" | (mu_correct & pi_correct))
  summary <- select(summary,
    seed, Version, Estimator, Mu, Pi, n_trial, n_aux,
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
    matches("^(truth|est|Var)_"),
    names_to = c("type", "parameter"), names_sep = "_"
  )
  summary <- pivot_wider(summary, names_from = type, values_from = value)
  summary <- mutate(summary,
    lower = est - 1.96 * sqrt(Var),
    upper = est + 1.96 * sqrt(Var)
  )
  summary <- group_by(summary, Version, Estimator, Mu, Pi, n_trial, n_aux, parameter)
  summary <- summarise(summary,
    bias = mean(est - truth),
    empirical_variance = var(est),
    estimated_variance = mean(Var),
    mse = mean((est - truth)^2),
    coverage = mean(truth >= lower & truth <= upper),
    .groups = "drop"
  )
  summary <- mutate(summary,
    Estimator = factor(Estimator, c("G-Formula", "IPW", "AIPW")),
    Version = factor(Version, c("Proposed", "Naive")),
    Pi = factor(Pi, c("-", "Incorrect Pi", "Correct Pi")),
    Mu = factor(Mu, c("-", "Incorrect Mu", "Correct Mu")),
    parameter = factor(
      parameter,
      c("eta0", "eta1", "rd", "rr"),
      c("eta(0)", "eta(1)", "RD", "RR")
    )
  )

  sample_size <- function(auxiliary, trial) {
    interaction(
      factor(auxiliary, levels = sort(unique(auxiliary))),
      factor(trial, levels = sort(unique(trial))),
      sep = ".", drop = TRUE
    )
  }
  formatted <- mutate(results,
    Estimator = factor(Estimator, c("G-Formula", "IPW", "AIPW")),
    Version = factor(Version, c("Proposed", "Naive")),
    Mu = if_else(
      Version == "Naive", "-",
      if_else(mu_correct, "Correct Mu", "Incorrect Mu")
    ),
    Pi = if_else(
      Version == "Naive", "-",
      if_else(pi_correct, "Correct Pi", "Incorrect Pi")
    ),
    consistent = Version == "Proposed" & (
      (Estimator == "G-Formula" & mu_correct) |
      (Estimator == "IPW" & pi_correct) |
      (Estimator == "AIPW" & (mu_correct | pi_correct))
    ),
    Sample_Size = sample_size(n_aux, n_trial)
  )
  formatted <- filter(formatted, Version == "Proposed" | (mu_correct & pi_correct))
  formatted <- mutate(formatted,
    Pi = factor(Pi, c("-", "Incorrect Pi", "Correct Pi")),
    Mu = factor(Mu, c("-", "Incorrect Mu", "Correct Mu"))
  )
  pi_labels <- c(
    "Correct Pi" = '"Correct "*pi[a]',
    "Incorrect Pi" = '"Incorrect "*pi[a]',
    "-" = '"-"'
  )
  mu_labels <- c(
    "Correct Mu" = '"Correct "*mu[a]',
    "Incorrect Mu" = '"Incorrect "*mu[a]',
    "-" = '"-"'
  )
  colors <- c("#FF6800", "#803E75", "#C10020", "#FFB300")
  shade <- distinct(formatted, Estimator, Version, Pi, Mu, consistent)
  shade <- filter(shade, !consistent)
  shade <- mutate(shade, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

  estimates <- ggplot(
    formatted,
    aes(Sample_Size, rdhat, color = Estimator, fill = Estimator)
  ) +
    geom_rect(
      data = shade,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#e5e4e2"
    ) +
    geom_boxplot(alpha = 0.5) +
    geom_hline(yintercept = mean(formatted$rd), linetype = "dashed") +
    facet_nested(
      Estimator ~ Version + Pi + Mu,
      labeller = labeller(
        Pi = as_labeller(pi_labels, label_parsed),
        Mu = as_labeller(mu_labels, label_parsed)
      )
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "Auxiliary and Trial Sample Sizes", y = expression(hat(RD))) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    guides(x = guide_axis_nested())

  variance_data <- filter(summary,
    n_trial == n_aux,
    empirical_variance > 0,
    estimated_variance > 0
  )
  variance <- ggplot(
    variance_data,
    aes(
      empirical_variance, estimated_variance,
      color = factor(n_trial), shape = parameter
    )
  ) +
    geom_rect(
      data = shade,
      aes(xmin = 0, xmax = xmax, ymin = 0, ymax = ymax),
      inherit.aes = FALSE, fill = "#e5e4e2"
    ) +
    geom_abline(linetype = "dashed") +
    geom_point(size = 3) +
    facet_nested(
      Estimator ~ Version + Pi + Mu,
      labeller = labeller(
        Pi = as_labeller(pi_labels, label_parsed),
        Mu = as_labeller(mu_labels, label_parsed)
      )
    ) +
    scale_x_continuous(
      transform = "log10",
      breaks = c(0.001, 0.01),
      labels = c("0.001", "0.01")
    ) +
    scale_y_continuous(
      transform = "log10",
      breaks = c(0.001, 0.01)
    ) +
    scale_color_manual(values = colors) +
    scale_shape_discrete(labels = function(x) parse(text = x)) +
    labs(
      x = "Empirical Variance", y = "Average Estimated Variance",
      color = "Trial and Auxiliary Sample Size", shape = "Parameter"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = unit(-5, "pt")
    )

  coverage_data <- mutate(summary, Sample_Size = sample_size(n_aux, n_trial))
  coverage <- ggplot(
    coverage_data,
    aes(Sample_Size, coverage, color = parameter, shape = parameter)
  ) +
    geom_rect(
      data = shade,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#e5e4e2"
    ) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    scale_color_manual(values = colors, labels = function(x) parse(text = x)) +
    scale_shape_discrete(labels = function(x) parse(text = x)) +
    facet_nested(
      Estimator ~ Version + Pi + Mu,
      labeller = labeller(
        Pi = as_labeller(pi_labels, label_parsed),
        Mu = as_labeller(mu_labels, label_parsed)
      )
    ) +
    labs(
      x = "Auxiliary and Trial Sample Sizes", y = "Empirical CI Coverage",
      color = "Estimand", shape = "Estimand"
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
    paste0("simulation1_", run_id, c("_estimates.png", "_variance.png", "_confidence.png"))
  )
  ggsave(files[1], estimates, width = 8, height = 6, dpi = dpi)
  ggsave(files[2], variance, width = 8, height = 6, dpi = dpi)
  ggsave(files[3], coverage, width = 8, height = 6, dpi = dpi)
  invisible(files)
}

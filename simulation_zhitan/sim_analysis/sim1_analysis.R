# Simulation 1 analysis ----------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(legendry)
library(grid)

sim1_colors <- c(
  "G-Formula" = "#FF6800", "IPW" = "#803E75",
  "AIPW" = "#C10020", "DML" = "#FFB300"
)

sim1_sample_size <- function(auxiliary, trial, separator = ".") {
  interaction(
    factor(auxiliary, levels = sort(unique(auxiliary))),
    factor(trial, levels = sort(unique(trial))),
    sep = separator, drop = TRUE
  )
}

sim1_size_labels <- function(labels, separator) {
  values <- strsplit(as.character(labels), separator, fixed = TRUE)
  as.expression(lapply(values, function(value) {
    bquote(n^aux == .(value[[1L]]) * "," ~~ n^trial == .(value[[2L]]))
  }))
}

sim1_summary <- function(results) {
  results %>%
    mutate(
      Pi = if_else(
        Version == "Naive", "-",
        if_else(pi_correct, "Correct Pi", "Incorrect Pi")
      ),
      Mu = if_else(
        Version == "Naive", "-",
        if_else(mu_correct, "Correct Mu", "Incorrect Mu")
      )
    ) %>%
    select(
      Version, Estimator, Mu, Pi, seed, n_trial, n_aux,
      eta_0, eta_1, rd, rr, etahat_0, etahat_1, rdhat, rrhat,
      cov_00, cov_11, var_rd, var_rr
    ) %>%
    rename(
      truth_eta0 = eta_0, truth_eta1 = eta_1,
      truth_rd = rd, truth_rr = rr,
      est_eta0 = etahat_0, est_eta1 = etahat_1,
      est_rd = rdhat, est_rr = rrhat,
      Var_eta0 = cov_00, Var_eta1 = cov_11,
      Var_rd = var_rd, Var_rr = var_rr
    ) %>%
    pivot_longer(
      cols = matches("^(truth|est|Var)_"),
      names_to = c("type", "parameter"), names_sep = "_"
    ) %>%
    pivot_wider(names_from = type, values_from = value) %>%
    mutate(
      lower = est - qnorm(0.975) * sqrt(Var),
      upper = est + qnorm(0.975) * sqrt(Var)
    ) %>%
    group_by(Version, Estimator, Mu, Pi, n_trial, n_aux, parameter) %>%
    summarise(
      empirical_variance = var(est),
      estimated_variance = mean(Var),
      bias = mean(est - truth),
      mse = mean((est - truth)^2),
      coverage = mean(truth >= lower & truth <= upper),
      .groups = "drop"
    ) %>%
    mutate(
      Estimator = factor(Estimator, c("G-Formula", "IPW", "AIPW", "DML")),
      Version = factor(Version, c("Proposed", "Naive")),
      Pi = factor(Pi, c("-", "Incorrect Pi", "Correct Pi")),
      Mu = factor(Mu, c("-", "Incorrect Mu", "Correct Mu")),
      parameter = factor(
        parameter, c("eta0", "eta1", "rd", "rr"),
        c("eta(0)", "eta(1)", "RD", "RR")
      )
    )
}

sim1_facet <- function(scales = "fixed") {
  facet_nested(
    Estimator ~ Version + Pi + Mu,
    scales = scales,
    labeller = labeller(
      Pi = as_labeller(c(
        "Correct Pi" = '"Correct "*pi[a]',
        "Incorrect Pi" = '"Incorrect "*pi[a]', "-" = '"-"'
      ), label_parsed),
      Mu = as_labeller(c(
        "Correct Mu" = '"Correct "*mu[a]',
        "Incorrect Mu" = '"Incorrect "*mu[a]', "-" = '"-"'
      ), label_parsed)
    )
  )
}

plot_sim1 <- function(
    results, figure_dir = "simulation/sim_figures/sim1", dpi = 600) {
  formatted <- results %>%
    mutate(
      Estimator = factor(Estimator, c("G-Formula", "IPW", "AIPW", "DML")),
      Version = factor(Version, c("Proposed", "Naive")),
      Pi = if_else(
        Version == "Naive", "-",
        if_else(pi_correct, "Correct Pi", "Incorrect Pi")
      ),
      Mu = if_else(
        Version == "Naive", "-",
        if_else(mu_correct, "Correct Mu", "Incorrect Mu")
      ),
      Pi = factor(Pi, c("-", "Incorrect Pi", "Correct Pi")),
      Mu = factor(Mu, c("-", "Incorrect Mu", "Correct Mu")),
      Sample_Size = sim1_sample_size(n_aux, n_trial),
      consistent = Estimator == "DML" | (Version == "Proposed" & (
        (Estimator == "G-Formula" & mu_correct) |
          (Estimator == "IPW" & pi_correct) |
          (Estimator == "AIPW" & (mu_correct | pi_correct))
      ))
    ) %>%
    filter(Version == "Proposed" | (mu_correct & pi_correct))

  summary <- sim1_summary(formatted)
  shade <- formatted %>%
    distinct(Estimator, Version, Pi, Mu, consistent) %>%
    filter(!consistent) %>%
    mutate(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  observed_dml <- formatted %>%
    filter(Estimator == "DML") %>%
    distinct(Estimator, Version, Pi, Mu)
  empty_dml <- formatted %>%
    distinct(Version, Pi, Mu) %>%
    mutate(Estimator = factor("DML", levels = levels(formatted$Estimator))) %>%
    anti_join(observed_dml, by = c("Estimator", "Version", "Pi", "Mu")) %>%
    mutate(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

  plot_theme <- theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    aspect.ratio = 1
  )
  estimates <- ggplot(
    formatted, aes(Sample_Size, rdhat, color = Estimator, fill = Estimator)
  ) +
    geom_rect(
      data = shade, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#C7C7C7"
    ) +
    geom_boxplot(alpha = 0.5, outlier.size = 0.35) +
    geom_hline(yintercept = mean(formatted$rd), linetype = "dashed") +
    geom_rect(
      data = empty_dml, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#666666"
    ) +
    sim1_facet("free_y") +
    scale_color_manual(values = sim1_colors) +
    scale_fill_manual(values = sim1_colors) +
    labs(x = "Auxiliary and Trial Sample Sizes", y = expression(hat(RD))) +
    plot_theme + theme(legend.position = "none") +
    guides(x = guide_axis_nested())

  variance_data <- summary %>%
    filter(empirical_variance > 0, estimated_variance > 0) %>%
    mutate(Sample_Size = sim1_sample_size(n_aux, n_trial, "_"))
  lower <- min(
    variance_data$empirical_variance, variance_data$estimated_variance
  ) / 2
  upper <- max(
    variance_data$empirical_variance, variance_data$estimated_variance
  ) * 2
  variance <- ggplot(
    variance_data,
    aes(empirical_variance, estimated_variance,
        color = Sample_Size, shape = parameter)
  ) +
    geom_rect(
      data = mutate(shade, xmin = lower, xmax = upper, ymin = lower, ymax = upper),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#C7C7C7"
    ) +
    geom_abline(linetype = "dashed") +
    geom_point(size = 3) +
    geom_rect(
      data = mutate(
        empty_dml, xmin = lower, xmax = upper, ymin = lower, ymax = upper
      ),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#666666"
    ) +
    sim1_facet() +
    scale_x_continuous(
      transform = "log10", breaks = c(0.001, 0.01),
      limits = c(lower, upper), expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      transform = "log10", breaks = c(0.001, 0.01),
      limits = c(lower, upper), expand = expansion(mult = 0)
    ) +
    scale_color_manual(
      values = unname(sim1_colors),
      labels = function(x) sim1_size_labels(x, "_")
    ) +
    scale_shape_discrete(labels = function(x) parse(text = x)) +
    labs(
      x = "Empirical Variance", y = "Average Estimated Variance",
      color = "Sample Size", shape = "Parameter"
    ) +
    plot_theme + theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.spacing.y = unit(-5, "pt")
    )

  coverage <- summary %>%
    mutate(Sample_Size = sim1_sample_size(n_aux, n_trial)) %>%
    ggplot(aes(Sample_Size, coverage, color = parameter, shape = parameter)) +
    geom_rect(
      data = shade, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#C7C7C7"
    ) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_rect(
      data = empty_dml, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#666666"
    ) +
    sim1_facet() +
    scale_color_manual(values = unname(sim1_colors)) +
    scale_shape_discrete(labels = function(x) parse(text = x)) +
    labs(
      x = "Auxiliary and Trial Sample Sizes", y = "Empirical CI Coverage",
      color = "Estimand", shape = "Estimand"
    ) +
    plot_theme + theme(legend.position = "bottom") +
    guides(x = guide_axis_nested())

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(
    figure_dir,
    c("sim1_estimates.png", "sim1_variance.png", "sim1_confidence.png")
  )
  ggsave(files[1], estimates, width = 9.5, height = 9.5, dpi = dpi)
  ggsave(files[2], variance, width = 9.5, height = 9.5, dpi = dpi)
  ggsave(files[3], coverage, width = 9.5, height = 9.5, dpi = dpi)
  invisible(files)
}

if (sys.nframe() == 0L) {
  sim1_files <- list.files(
    "simulation/sim_data/sim1", "^sd[0-9]+[.]csv$", full.names = TRUE
  )
  sim1_results <- bind_rows(lapply(sim1_files, read.csv))
  plot_sim1(sim1_results)
}

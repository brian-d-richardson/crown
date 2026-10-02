# Simulation 2 analysis ----------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(legendry)
library(grid)

sim2_colors <- c(
  "G-Formula" = "#FF6800", "IPW" = "#803E75",
  "AIPW" = "#C10020", "DML" = "#FFB300"
)

sim2_sample_size <- function(auxiliary, trial, separator = ".") {
  interaction(
    factor(auxiliary, levels = sort(unique(auxiliary))),
    factor(trial, levels = sort(unique(trial))),
    sep = separator, drop = TRUE
  )
}

sim2_size_labels <- function(labels, separator) {
  values <- strsplit(as.character(labels), separator, fixed = TRUE)
  as.expression(lapply(values, function(value) {
    bquote(n^aux == .(value[[1L]]) * "," ~~ n^trial == .(value[[2L]]))
  }))
}

sim2_summary <- function(results) {
  results %>%
    select(
      Estimator_Panel, seed, n_trial, n_aux,
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
    group_by(Estimator_Panel, n_trial, n_aux, parameter) %>%
    summarise(
      empirical_variance = var(est),
      estimated_variance = mean(Var),
      bias = mean(est - truth),
      mse = mean((est - truth)^2),
      coverage = mean(truth >= lower & truth <= upper),
      .groups = "drop"
    ) %>%
    mutate(
      Estimator_Panel = factor(
        Estimator_Panel,
        c("Proposed G-Formula", "Proposed IPW", "Proposed AIPW", "Proposed DML")
      ),
      parameter = factor(
        parameter, c("eta0", "eta1", "rd", "rr"),
        c("eta(0)", "eta(1)", "RD", "RR")
      )
    )
}

plot_sim2 <- function(
    results, figure_dir = "simulation/sim_figures/sim2", dpi = 600) {
  labels <- c(
    "G-Formula" = "Proposed G-Formula", "IPW" = "Proposed IPW",
    "AIPW" = "Proposed AIPW", "DML" = "Proposed DML"
  )
  formatted <- results %>%
    filter(Version == "Proposed", Estimator %in% names(labels)) %>%
    mutate(
      Estimator = factor(Estimator, names(labels)),
      Estimator_Panel = factor(labels[as.character(Estimator)], labels),
      Sample_Size = sim2_sample_size(n_aux, n_trial)
    )
  summary <- sim2_summary(formatted)
  plot_theme <- theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

  estimates <- ggplot(
    formatted, aes(Sample_Size, rdhat, color = Estimator, fill = Estimator)
  ) +
    geom_boxplot(alpha = 0.5, outlier.size = 0.35) +
    geom_hline(yintercept = mean(formatted$rd), linetype = "dashed") +
    facet_nested(. ~ Estimator_Panel) +
    scale_color_manual(values = sim2_colors) +
    scale_fill_manual(values = sim2_colors) +
    labs(x = "Auxiliary and Trial Sample Sizes", y = expression(hat(RD))) +
    plot_theme + theme(legend.position = "none") +
    guides(x = guide_axis_nested())

  variance <- summary %>%
    filter(empirical_variance > 0, estimated_variance > 0) %>%
    mutate(Sample_Size = sim2_sample_size(n_aux, n_trial, "_")) %>%
    ggplot(aes(
      empirical_variance, estimated_variance,
      color = Sample_Size, shape = parameter
    )) +
    geom_abline(linetype = "dashed") +
    geom_point(size = 3) +
    facet_nested(. ~ Estimator_Panel) +
    scale_x_continuous(transform = "log10", breaks = c(0.001, 0.01)) +
    scale_y_continuous(transform = "log10", breaks = c(0.001, 0.01)) +
    scale_color_manual(
      values = unname(sim2_colors),
      labels = function(x) sim2_size_labels(x, "_")
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
    mutate(Sample_Size = sim2_sample_size(n_aux, n_trial)) %>%
    ggplot(aes(Sample_Size, coverage, color = parameter, shape = parameter)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    facet_nested(. ~ Estimator_Panel) +
    scale_color_manual(values = unname(sim2_colors)) +
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
    c("sim2_estimates.png", "sim2_variance.png", "sim2_confidence.png")
  )
  ggsave(files[1], estimates, width = 12, height = 3.7, dpi = dpi)
  ggsave(files[2], variance, width = 12, height = 4.8, dpi = dpi)
  ggsave(files[3], coverage, width = 12, height = 4.8, dpi = dpi)
  invisible(files)
}

if (sys.nframe() == 0L) {
  sim2_files <- list.files(
    "simulation/sim_data/sim2", "^sd[0-9]+[.]csv$", full.names = TRUE
  )
  sim2_results <- bind_rows(lapply(sim2_files, read.csv))
  plot_sim2(sim2_results)
}

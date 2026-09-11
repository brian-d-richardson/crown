#!/usr/bin/env Rscript

# Simulation 2: all supported estimator combinations

project_dir <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", ".."))
setwd(project_dir)

# Packages

library(dplyr)
library(tidyr)
library(ggplot2)
library(legendry)
library(grid)
library(crown)

# Functions

source("simulation/sim_functions/sim2.R")
source("simulation/sim_analysis/sim2_analysis.R")

# Settings

sample_sizes <- expand.grid(
  n_trial = c(500L, 5000L),
  n_auxiliary = c(500L, 5000L)
)
mc_reps <- 20L

run_id <- "all_estimators_mc20"
data_dir <- file.path(project_dir, "simulation", "sim_data", "sim2")
figure_dir <- file.path(project_dir, "simulation", "sim_figures", "sim2")

# Run

simulation <- run_sim2(
  sample_sizes = sample_sizes,
  run_id = run_id,
  out_dir = data_dir,
  mc_reps = mc_reps
)

figures <- plot_sim2(
  simulation$results,
  run_id,
  figure_dir
)

cat(
  "\nCompleted", mc_reps, "Monte Carlo replicates for each of",
  nrow(sample_sizes), "sample-size combinations\n"
)
cat("Results:", simulation$result_file, "\n")
cat("Figures:\n")
print(figures)

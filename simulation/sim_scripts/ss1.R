#!/usr/bin/env Rscript

# Simulation 1: parametric G-formula, IPW, and AIPW

project_dir <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", ".."))
setwd(project_dir)

# Packages

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(legendry)
library(grid)
library(crown)

# Functions

source("simulation/sim_functions/sim1.R")
source("simulation/sim_analysis/sim1_analysis.R")

# Settings

sample_sizes <- expand.grid(
  n_trial = c(500L, 5000L),
  n_auxiliary = c(500L, 5000L)
)
mc_reps <- 20L
n_clusters <- 20L
p_response <- 0.5
p_censoring <- 0.3
run_id <- "parametric_mc20"
data_dir <- file.path(project_dir, "simulation", "sim_data", "sim1")
figure_dir <- file.path(project_dir, "simulation", "sim_figures", "sim1")

# Run

simulation <- run_sim1(
  sample_sizes = sample_sizes,
  run_id = run_id,
  out_dir = data_dir,
  n_clusters = n_clusters,
  p_response = p_response,
  p_censoring = p_censoring,
  mc_reps = mc_reps
)

figures <- plot_sim1(
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

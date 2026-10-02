# Simulation 2 -------------------------------------------------------------

# Setup

project_dir <- normalizePath(
  file.path(dirname(sys.frame(1)$ofile), "..", "..")
)
setwd(project_dir)

devtools::load_all()
source("simulation/sim_functions/sim2.R")
source("simulation/sim_analysis/sim2_analysis.R")

# Settings

# Local example: 50 Monte Carlo replicates.
# HPC: 10 independent jobs x 1,000 replicates = 10,000 total replicates.
# The HPC shell script should assign a different seed and output file to each job.
mc_reps <- 50

sample_sizes <- expand.grid(
  n_trial = c(500, 5000),
  n_auxiliary = c(500, 5000)
)

# Run

simulation <- run_sim2(
  sample_sizes,
  run_id = "sd_local",
  out_dir = "simulation/sim_data/sim2",
  mc_reps = mc_reps,
  base_seed = 91000000,
  arguments = list(nthread = 1)
)

figures <- plot_sim2(
  simulation$results,
  figure_dir = "simulation/sim_figures/sim2/local"
)

###############################################################################
###############################################################################

# PopART Simulation 2 Script

# Brian Richardson

# 2026-08-24

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

# clear workspace
rm(list = ls())

# load necessary packages
library(dplyr)
library(tidyr)
library(devtools)
library(xgboost)
library(pbapply)

# check if running on cluster (env variable is non-empty)
env_task <- Sys.getenv("SLURM_ARRAY_TASK_ID")
on.cluster <- env_task != ""
if (on.cluster) {
  cluster.id <- as.numeric(env_task)
  setwd(dirname(getwd()))
} else {
  cluster.id <- 0
  setwd("C:/Users/brich/OneDrive - University of North Carolina at Chapel Hill/Desktop/CIRL/PopART/crown")
}

# load crown and simulation functions
load_all()
source("simulation/sim_functions/sim2.R")

# simulation parameters ---------------------------------------------------

## baseline seed (specific to cluster)
base.seed <- 10^6 * as.integer(cluster.id)

## fixed parameters
m <- 20
p_resp <- 0.5
p_cens <- 0.3

## number of simulation replicates
n.rep <- 1

## simulation inputs
sim.in <- expand.grid(
  n_trial = c(500, 5000),
  n_aux = c(500, 5000),
  K = c(5, 10),
  sim.id = 1:n.rep + base.seed)

## test run one simulation
if (FALSE) {
  sim2_fun(
    m = m,
    n_trial = sim.in$n_trial[1],
    n_aux = sim.in$n_aux[1],
    K = sim.in$K[1],
    p_resp = p_resp,
    p_cens = p_cens,
    seed = sim.in$sim.id[1])
}

# run simulations ---------------------------------------------------------

## run simulations
sim.out <- pblapply(
  X = seq_len(nrow(sim.in)),
  FUN = function(ii) {
    out <- tryCatch(
      {
        sim2_fun(
          m = m,
          n_trial = sim.in$n_trial[ii],
          n_aux = sim.in$n_aux[ii],
          K = sim.in$K[ii],
          p_resp = p_resp,
          p_cens = p_cens,
          seed = sim.in$sim.id[ii])
      },
      error = function(e) {
        message(sprintf("Error in row %d (sim.id = %s): %s",
                        ii, sim.in$sim.id[ii], conditionMessage(e)))
        NULL
      }
    )
    out
  }) %>%
  bind_rows()

## save results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim2/sd",
                 as.integer(cluster.id), ".csv"))

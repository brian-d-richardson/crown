###############################################################################
###############################################################################

# PopART Simulation 1 Script

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
source("simulation/sim_functions/sim1.R")

# simulation parameters ---------------------------------------------------

## baseline seed (specific to cluster)
base.seed <- 10^6 * as.integer(cluster.id)

## fixed parameters
m <- 20
p_resp <- 0.5
p_cens <- 0.3

## number of simulation replicates
n.rep <- 2

## simulation inputs
sim.in <- expand.grid(
  n_trial = c(500, 5000),
  n_aux = c(500, 5000),
  mu_correct = c(T, F),
  pi_correct = c(T, F),
  sim.id = 1:n.rep + base.seed)

## test run one simulation
if (FALSE) {
  sim1_fun(
    m = m,
    n_trial = sim.in$n_trial[1],
    n_aux = sim.in$n_aux[1],
    mu_correct = sim.in$mu_correct[1],
    pi_correct = sim.in$pi_correct[1],
    p_resp = p_resp,
    p_cens = p_cens,
    seed = sim.in$sim.id[1])
}

# run simulations ---------------------------------------------------------

## run simulations
sim.out <- pblapply(
  X = seq_len(nrow(sim.in)),
  FUN = function(ii) {

    sim1_fun(
      m = m,
      n_trial = sim.in$n_trial[ii],
      n_aux = sim.in$n_aux[ii],
      mu_correct = sim.in$mu_correct[ii],
      pi_correct = sim.in$pi_correct[ii],
      p_resp = p_resp,
      p_cens = p_cens,
      seed = sim.in$sim.id[ii])

  }) %>%
  bind_rows()

## save results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim1/sd",
                 as.integer(cluster.id), ".csv"))

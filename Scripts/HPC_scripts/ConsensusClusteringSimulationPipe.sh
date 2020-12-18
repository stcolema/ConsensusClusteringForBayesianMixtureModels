#!/bin/bash

First create the data
Rscript /home/sdc56/rds/hpc-work/Scripts/ConsensusInference/R/scenarios.R



# Call MDI and Mclust
### MDI call ###
### Scenario call script ###
# NOTE: this is probably infeasible and various jobs will die/fail
bash ./SLURM/all_scenarios.sh

# CC model eval
bash ./R/CallCIModelEval.R

# Check model convergence for Bayesian
## Convergence within Sim
## Convergence across Sim

# Bayes model evaluation

# Model evaluation
## comparing all models


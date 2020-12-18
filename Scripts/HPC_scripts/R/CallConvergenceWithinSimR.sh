#!/bin/bash

scenarios=(
"simple_2d"
"no_structure"
"base_case"
"large_standard_deviation_3"
"large_standard_deviation_5"
"irrelevant_features_10"
"irrelevant_features_20"
"irrelevant_features_100"
"varying_proportions"
"large_n_small_p_base"
"large_n_small_p_large_k"
"large_n_small_p_large_k_small_dm"
"small_n_large_p_base"
"small_n_large_p_small_dm"
"varying_proportions_small_dm"
"large_n_large_p"
"large_n_large_p_small_dm"
);


# scn="no_structure";
for scn in "${scenarios[@]}"; do
  for i in {1..100}; do 
    Rscript ConvergenceWithinSim.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${i}/Bayesian/ -s ~/rds/hpc-work/Analysis/Simulations/Convergence/${scn}/Within_chain/ -l 1000001 &
  done;
  sleep 240
done;

#!/bin/bash

# perform model analysis of separate CI runs within a scenario
# simulations 09/05/2020

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
# "large_n_small_p_base"
# "large_n_small_p_large_k"
# "large_n_small_p_large_k_small_dm"
"small_n_large_p_base"
"small_n_large_p_small_dm" 
"varying_proportions_small_dm"
#"large_n_large_p"
#"large_n_large_p_small_dm"
);


# scn="base_case";
# scn="large_n_small_p_base"

for scn in "${scenarios[@]}"; do 
  for s in {1..100}; do 
    Rscript CIModelEval.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${s}/Consensus/ --truth_dir ~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/${scn}/ --sim_num ${s} --save_dir ~/rds/hpc-work/Analysis/Simulations/Model_performance/${scn}/Consensus/ --scn ${scn} &
  done;
  sleep 250;
done;

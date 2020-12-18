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
);


# scn="no_structure";
# nice_scn="No Structure";


for scn in "${scenarios[@]}"; do
  Rscript ConvergenceAcrossSim.R -d  ~/rds/hpc-work/Analysis/Simulations/Convergence/${scn}/Within_chain/ -s ~/rds/hpc-work/Analysis/Simulations/Convergence/${scn}/Across_chain/ --scn ${scn} -f ${scn}ConvergenceAcrossChains.png &
done; 

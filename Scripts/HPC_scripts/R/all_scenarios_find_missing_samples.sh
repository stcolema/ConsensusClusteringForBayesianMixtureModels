#!/usr/env/bin/bash

scenarios=("simple_2d"
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
)

n_scn=${#scenarios[@]}

for (( i=0; i<${n_scn}; i++ ));
do
  for (( j=0; j<${n_sim}; j++ ));
  do
    echo ${scenarios[$i]};
    Rscript ./FindMissingSamples.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${j}/Consensus/ > missing${scn}Sim${j}.txt
  done;
done;

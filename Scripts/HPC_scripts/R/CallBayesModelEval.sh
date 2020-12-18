#!/bin/bash


scn="no_structure"; 


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

for scn in "${scenarios[@]}"; do
for s in {1..100}; do 
  Rscript BayesModelEval.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${s}/Bayesian/ --truth_dir ~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/${scn}/ --sim_num ${s} --save_dir ~/rds/hpc-work/Analysis/Simulations/Model_performance/${scn}/Bayesian/ --geweke_dir ~/rds/hpc-work/Analysis/Simulations/Convergence/${scn}/Within_chain/ --scn ${scn} & 
done;
sleep 240;
done;


#!/usr/env/bin/bash

scn="simple_2d"; 
n_sim=100;

for i in {1..${n_sim}}; do 
  Rscript ~/hpc_stuff/Consensus_clustering/scripts/R_scripts/Simulations/CIModelEval.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${i}/Consensus/ --sim_num ${i} --save_dir ~/rds/hpc-work/Analysis/Simulations/Model_performance/${scn}/Consensus/ --truth_dir ~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/${scn}/ --cols "1"; 
done


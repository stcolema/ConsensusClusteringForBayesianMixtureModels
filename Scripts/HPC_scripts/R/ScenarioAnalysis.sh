#!/bin/bash

# Bash script to implement full analysis pipeline of simulations to compare
# MLE, Bayesian inference and Consensus inference of mixture models.

# Base case
scn="base_case";


# Receive parameters from command line
while getopts s: option; do
    case "${option}" in
        s) scn=${OPTARG};;
    esac
done


# Bayesian inference:
# Check convergence both within and across simulations 
echo "Check convergence both within and across chains for each simulation." 
for i in {1..100}; do
  Rscript ConvergenceWithinSim.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${i}/Bayesian/ -s ~/rds/hpc-work/Analysis/Simulations/Convergence/${scn}/Within_chain/ -l 1000001 &
  pids[${i}]=$!
done

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done


echo "Plot Gelman-Rubin values across all simulations to enable comparison." 
Rscript ConvergenceAcrossSim.R -d  ~/rds/hpc-work/Analysis/Simulations/Convergence/${scn}/Within_chain/ -s ~/rds/hpc-work/Analysis/Simulations/Convergence/${scn}/Across_chain/ --scn ${scn} -f ${scn}ConvergenceAcrossChains.png 



# Evaluate model performance
echo "Evaluate model performance compared to the true labelling used to generate the data."

# 1. Consensus inference
echo "Consensus inference"
for s in {1..100}; do
  Rscript CIModelEval.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${s}/Consensus/ --truth_dir ~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/${scn}/ --sim_num ${s} --save_dir ~/rds/hpc-work/Analysis/Simulations/Model_performance/${scn}/Consensus/ --scn ${scn} &
 pids[${s}]=$!
done

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done


# 2. Bayesian inference
echo "Bayesian inference"
for s in {1..100}; do
  Rscript BayesModelEval.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${s}/Bayesian/ --truth_dir ~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/${scn}/ --sim_num ${s} --save_dir ~/rds/hpc-work/Analysis/Simulations/Model_performance/${scn}/Bayesian/ --geweke_dir ~/rds/hpc-work/Analysis/Simulations/Convergence/${scn}/Within_chain/ --scn ${scn} &
 pids[${s}]=$!
done

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done;


# 3. MLE (Mclust)
echo "MLE / Mclust".
for sim_num in {1..100}; do 
  Rscript FreqModelEval.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${sim_num}/Frequentist/ -s ~/rds/hpc-work/Analysis/Simulations/Model_performance/${scn}/Frequentist/ --sim_num ${sim_num} --truth_dir ~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/${scn}/ &
 pids[${sim_num}]=$!
done

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done;


# Compare all models!
echo "Compare all methods creating plots to enable interpretation."
Rscript ComparingAllModels.R -d ~/rds/hpc-work/Analysis/Simulations/Model_performance/ --scn "${scn}" --save_dir ~/rds/hpc-work/Analysis/Simulations/Model_performance/${scn}/ --truth_dir ~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/${scn}/

echo "Done."

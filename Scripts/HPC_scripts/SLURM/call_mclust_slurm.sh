#!/usr/env/bin/bash

# example: bash call_mclust_slurm.sh -o "small_n_large_p_base" -i 1 -k 2 -g 25;


# Directory to save output to
save_dir="/home/sdc56/rds/hpc-work/MDI_output/Simulations/Single_dataset/";

# Random seed to use
seed=1;

# Directory to read data from
curr_dir="/home/sdc56/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/";

# Scenario to run
scn="base_case";

# Number of clusters allowed in MDI
n_clust=2;
n_clust_upper=50;

# Simulation number
sim_num=1;

# Inference type; one of Consensus (default), Bayesian or Frequentist
model="Frequentist";

# Model name for Mclust
model_name="NULL";

# Receive parameters from command line
while getopts o:i:k:g:m: option; do
    case "${option}" in
        o) scn=${OPTARG};;
        i) sim_num=${OPTARG};;
        k) n_clust=${OPTARG};;
        g) n_clust_upper=${OPTARG};;
        m) model_name=${OPTARG};;
    esac
done

# Output and input homes
input_data="${curr_dir}${scn}/dataset_${sim_num}.csv";
output_dir="${save_dir}${scn}/simulation_${sim_num}/${model}/";

# The output file name
output="${output_dir}${model}ModType${model_name}.Rds";

# File to save timing to
time_file="./Time/${scn}/${model}/time_scn_${scn}_sim_${sim_num}_model_${model}$M{model_name}.txt";

# Application to run
application="Rscript ~/rds/hpc-work/Scripts/ConsensusInference/R/call_mclust.R";

# Arguments to pass to application
options=" -d ${input_data} -k ${n_clust} -g ${n_clust_upper} -m ${model_name} -f ${output}";

# Combine these and add timing call
CMD="{ time $application $options; } 2> ${time_file} ";

# Print and call command
# echo ${CMD};
eval ${CMD};


#!/usr/bin/bash

source /home/sdc56/rds/hpc-work/Scripts/ConsensusInference/SLURM/CIFindTime.sh

x1="$(calculate_sbatch_time 1000001 16)";
echo  "$x1";

x2="$(calculate_sbatch_time 1000001 175)";
echo "${x2}";

x3="$(calculate_sbatch_time 1000001 1)";
echo $x3

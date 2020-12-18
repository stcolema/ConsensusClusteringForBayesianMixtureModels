#!/bin/bash


function calculate_sbatch_time {

num_iter=$1;
iter_5k=$2;


let num_sec=${num_iter}*${iter_5k}/5000;
let num_min=${num_sec}/60;
let num_hours=${num_min}/60;

if [ $num_hours -ge 12 ]
then
  num_hours=11;
  num_min=719;
  num_sec=43199;
fi

# Adjust
let num_min=${num_min}-60*${num_hours};
let num_sec=${num_sec}-60*${num_min}-3600*${num_hours};

# The string length of each of these
hour_str_length=${#num_hours};
min_str_length=${#num_min};
sec_str_length=${#num_sec};

# If any length == 1 add some 0
hour_str="${num_hours}";

if [ ${hour_str_length} -lt 2 ]; then
  hour_str="0${hour_str}";
fi

min_str="${num_min}";
if [ ${min_str_length} -lt 2 ]; then
  min_str="0${min_str}";
fi

sec_str="${num_sec}";
if [ ${sec_str_length} -lt 2 ]; then
  sec_str="0${sec_str}";
fi

# File to save timing to
time_file="/home/sdc56/rds/hpc-work/Scripts/ConsensusInference/SLURM/Time/${scn}/${model}/${scn}_sim_${sim_num}_model_${model}_N_${num_iter}_seed_${seed}_K_${n_clust}.txt"

time_str="${hour_str}:${min_str}:${sec_str}";

echo "$time_str";
}

#time_1="$(calculate_sbatch_time 10001 15)";
#time_2="$(calculate_sbatch_time 1000001 15)";
#
#echo "$time_1";
#echo "$time_2"

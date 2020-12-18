#!/bin/bash

source /home/sdc56/rds/hpc-work/Scripts/ConsensusInference/SLURM/CIFindTime.sh

# Scenario run
scn="simple_2d";
time_for_5000=5;
n_sim=100;
n_c=100;
n_b=10;
n_f=1;
n_iter_b=1000001;
thin_b=1000;
k=50;

# Receive parameters from command line
while getopts s:t:n:c:b:f:k: option; do
    case "${option}" in
        s) scn=${OPTARG};;
        t) time_for_5000=${OPTARG};;
        n) n_sim=${OPTARG};;
        c) n_c=${OPTARG};;
        b) n_b=${OPTARG};;
        f) n_f=${OPTARG};;
        k) k=${OPTARG};;
    esac
done

# Directory to save output to
save_dir="/home/sdc56/rds/hpc-work/MDI_output/Simulations/Single_dataset/";

large_k_scn=("large_n_small_p_base"
"large_n_small_p_large_k"
"large_n_small_p_large_k_small_dm"
)

check_sim=75;
if [[ ${check_sim} -gt ${n_sim} ]]; then
  check_sim=${n_sim};
fi

# Use 10001 as 0-10,000 used
for n in 10 100 1000 10001; do
  let thin=${n}/10; 

  # iterate over random seed / chain number
  for ((s = 1; s <= ${n_c}; s++)); do
    if [ ${n} -le 1000 ]
    then 
      for ((i = 1; i <= ${n_sim}; i++)); do
      
      output_dir="${save_dir}${scn}/simulation_${i}/Consensus/";
  
      # The output file name
      file_check="${output_dir}ConsensusN${n}T${thin}Seed${s}.csv";
      file_check2="${output_dir}ConsensusN${n}T${thin}Seed${s}K${k}.csv";

      # If in a large K case the file name of importane is different
      if [[ ${k} -ne 50 ]]; then
        # The output file name
        file_check="${output_dir}ConsensusN${n}T${thin}Seed${s}K${k}.csv";
      fi
      

      if [[ -f "$file_check" ]] || [[ -f "$file_check2" ]]; then
        echo "$file_check exists, skipping run."
        continue
      fi


        bash run_mason.sh -s ${s} -i ${i} -n ${n} -m Consensus -t ${thin} -o ${scn} -l ${time_for_5000} -h 0 -k ${k}; 
      done;
    else



      output_dir="${save_dir}${scn}/simulation_${check_sim}/Consensus/";

      # The output file name
      file_check="${output_dir}ConsensusN${n}T${thin}Seed${s}.csv";
      file_check2="${output_dir}ConsensusN${n}T${thin}Seed${s}K${k}.csv";

      # If in a large K case the file name of importane is different
      if [[ " ${large_k_scn[@]} " =~ " ${scn} " ]]; then
        # The output file name
        file_check="${output_dir}ConsensusN${n}T${thin}Seed${s}K${k}.csv";
      fi

      ###if [[ -f "$file_check" ]]; then
      if [[ -f "$file_check" ]] || [[ -f "$file_check2" ]]; then
        echo "$file_check exists, skipping run."
        continue;
        ###exit 0
      fi

      time_str="$(calculate_sbatch_time ${n} ${time_for_5000})";
      job_name="${scn}ConsensusN${n}";
      
      sbatch --job-name=${job_name} --array=1-${n_sim} --time=${time_str} run_mason_array.sh -s ${s} -n ${n} -m Consensus -t ${thin} -o ${scn} -l ${time_for_5000} -k ${k};
    fi;
  done & 
done;

# Bayesian inference
n=${n_iter_b};
thin=${thin_b};

# as above, iterate over simulation (redundant) and chain
for ((s = 1; s <= ${n_b}; s++)); do
  ##bash run_mason.sh -s ${s} -i ${i} -n ${n} -m Bayesian -t ${thin} -o ${scn} -l ${time_for_5000} -h 0;

  output_dir="${save_dir}${scn}/simulation_${check_sim}/Bayesian/";

  # The output file name
  file_check="${output_dir}BayesianN${n}T${thin}Seed${s}.csv";
  file_check2="${output_dir}BayesianN${n}T${thin}Seed${s}K${k}.csv";

  if [[ -f "$file_check" ]] || [[ -f "${file_check2}" ]] ; then
    echo "$file_check exists, skipping run."
    continue;
  fi

  time_str="$(calculate_sbatch_time ${n} ${time_for_5000})";
  job_name="${scn}BayesianN${n}";
 
  echo "${time_str}";


  sbatch --job-name=${job_name} --time="${time_str}" --array=1-${n_sim} run_mason_array.sh -s ${s} -n ${n} -m Bayesian -t ${thin} -o ${scn} -l ${time_for_5000} -k ${k};
done;

# Frequentist inference
for ((i = 0; i <= ${n_sim}; i++)); do 
  bash call_mclust_slurm.sh -o ${scn} -i  ${i}; 
done;

# Evaluate Consensus inference model performance
#for i in {1..100}; do
  ###Rscript ~/hpc_stuff/Consensus_clustering/scripts/R_scripts/Simulations/CIModelEval.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${i}/Consensus/ --sim_num ${i} --save_dir ~/rds/hpc-work/Analysis/Simulations/Model_performance/${scn}/Consensus/ --truth_dir ~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/${scn}/ --cols "1";
# done


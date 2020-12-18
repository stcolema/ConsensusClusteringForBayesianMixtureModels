
scn="base_case";
scn="simple_2d";
scn="no_structure";
scn="large_standard_deviation_3";
scn="large_standard_deviation_5";
scn="irrelevant_features_10";
scn="irrelevant_features_20";
scn="irrelevant_features_100";
scn="varying_proportions";
scn="large_n_small_p_base";
scn="large_n_small_p_large_k";
scn="large_n_small_p_large_k_small_dm"
scn="small_n_large_p_base";
scn="small_n_large_p_small_dm";

scenarios=(
# "simple_2d"
# "no_structure"
# "base_case"
# "large_standard_deviation_3"
# "large_standard_deviation_5"
# "irrelevant_features_10"
# "irrelevant_features_20"
"irrelevant_features_100"
# "varying_proportions"
# "large_n_small_p_base"
# "large_n_small_p_large_k"
# "large_n_small_p_large_k_small_dm"
# "small_n_large_p_base"
# "small_n_large_p_small_dm"
);



inf="Consensus";
n=10001;
t=1000;
k=50;

sim_l=1;
sim_u=101;

seed_l=1;
seed_u=101;

for scn in "${scenarios[@]}"; do
  for (( i=${sim_l}; i < ${sim_u}; i++)); do 
    # echo "Sim ${i}"; 
    Rscript FindMissingSamples.R -d ~/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${i}/${inf}/ -s "${scn}Missing${inf}Sim${i}.csv" -i ${n} -t ${t} -n 100 --inference_type ${inf} -k ${k} &
    pids[${i}]=$!;
  done;

  # wait for all pids
  for pid in ${pids[*]}; do
    wait $pid;
  done;
done;





scn="large_n_small_p_base";
n=10001;
t=1000;
k=50;

sim_l=1;
sim_u=101;

seed_l=111111111111111111111111;
seed_u=101;

for (( i = ${seed_l}; i < ${seed_u}; i++ )); do
  for (( s = ${sim_l}; s < ${sim_u}; s++)); do
    if [ ! -f "/home/sdc56/rds/hpc-work/MDI_output/Simulations/Single_dataset/${scn}/simulation_${s}/Consensus/ConsensusN${n}T${t}Seed${i}K${k}.csv" ]; then 
      bash run_mason.sh -o ${scn} -l 15 -k ${k} -t ${t} -n ${n} -s ${i} -i ${s} -m "Consensus" & 
    fi;
  done; 
done



files_to_drop=(ConsensusN10001T1000Seed25K50.csv 
ConsensusN10001T1000Seed26K50.csv
ConsensusN10001T1000Seed27K50.csv
ConsensusN10001T1000Seed28K50.csv
ConsensusN10001T1000Seed30K50.csv
ConsensusN10001T1000Seed34K50.csv
ConsensusN10001T1000Seed35K50.csv
ConsensusN10001T1000Seed36K50.csv
ConsensusN10001T1000Seed42K50.csv
ConsensusN10001T1000Seed50K50.csv
ConsensusN10001T1000Seed53K50.csv
ConsensusN10001T1000Seed55K50.csv
ConsensusN10001T1000Seed56K50.csv
);

for f in "${files_to_drop[@]}"; do 
  rm "/home/sdc56/rds/hpc-work/MDI_output/Simulations/Single_dataset/small_n_large_p_base/simulation_1/Consensus/${f}"; 
done;

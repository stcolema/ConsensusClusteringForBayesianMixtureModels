#!/usr/env/bin/bash

scenarios=("simple_2d"
"no_structure"
"base_case"
"large_standard_deviation_3"
"large_standard_deviation_5"
"irrelevant_features_10"
"irrelevant_features_20"
"irrelevant_features_100"
"varying_proportions"
"large_n_small_p_base"
"large_n_small_p_large_k"
"large_n_small_p_large_k_small_dm"
"small_n_large_p_base"
"small_n_large_p_small_dm"
"varing_proportions_small_dm"
"large_n_large_p"
"large_n_large_p_small_dm"
)

#The time for 5000 iterations in each scenario (ordered as above)
times=(7
7
12
12
12
18
20
56
12
175 # 1690
175 # 1690
175 # 1690
75
75
12
175
175
)

# The number of clusters for each scenario
k=(50
50
50
50
50
50
50
50
50
50
50
50
50
50
50
250
250
)


n_scn=${#scenarios[@]}

for (( i=0; i<${n_scn}; i++ ));
do
  echo ${scenarios[$i]};
  run_scenario.sh -s ${scenarios[$i]} -t ${times[$i]} -k ${k[$i]} &
done

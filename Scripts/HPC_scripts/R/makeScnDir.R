
scenarios=c("simple_2d",
  "no_structure",
  "base_case",
"large_standard_deviation_3",
"large_standard_deviation_5",
"irrelevant_features_10",
"irrelevant_features_20",
"irrelevant_features_100",
"varying_proportions",
"large_n_small_p_base",
"large_n_small_p_large_k",
"large_n_small_p_large_k_small_dm",
"small_n_large_p_base",
"small_n_large_p_small_dm"
);

inf=c("Bayesian", "Consensus", "Frequentist")


main_dir="~/rds/hpc-work/Analysis/Simulations/Model_performance/"
for(scn in scenarios){
  scn_dir=paste0(main_dir, scn)
  dir.create(scn_dir, showWarnings = FALSE)
  for(i in inf){
    inf_dir= paste0(scn_dir, "/", i)
    dir.create(inf_dir, showWarnings = FALSE)
  }
}


types=c("Across_chain", "Within_chain");
main_dir="~/rds/hpc-work/Analysis/Simulations/Convergence/";
for(scn in scenarios){
  scn_dir=paste0(main_dir, scn)
  dir.create(scn_dir, showWarnings = FALSE)
  for(t in types){
    inf_dir= paste0(scn_dir, "/", t)
    dir.create(inf_dir, showWarnings = FALSE)
  }
}




scn="small_n_large_p_small_dm";


# Receive parameters from command line
while getopts s: option; do
    case "${option}" in
        s) scn=${OPTARG};;
    esac
done


Rscript ComparingAllModels.R -d ~/rds/hpc-work/Analysis/Simulations/Model_performance/ --scn "${scn}" --save_dir ~/rds/hpc-work/Analysis/Simulations/Model_performance/${scn}/ --truth_dir ~/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/${scn}/

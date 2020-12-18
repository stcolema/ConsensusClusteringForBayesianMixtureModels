#!/bin/bash
#!
#! Example call: sbatch --job-name=${job_name} --array=1-${n_sim} --time=${time_str} run_mason_array.sh -s ${s} -n ${n} -m Consensus -t ${thin} -o ${scn} -l ${time_for_5000} -k ${k}

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
##SBATCH -J run_simulation_model 

#! Which project should be charged:
##SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -A CWALLACE-SL3-CPU


#! How many whole nodes should be allocated?
#SBATCH --nodes=1

#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1

#! How much wallclock time will be required?
##SBATCH --time=00:00:15

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem": 
#SBATCH -p skylake-himem

#! sbatch directives end here (put any additional directives above this line)

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numtasks=$SLURM_NTASKS
numnodes=$SLURM_JOB_NUM_NODES
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

## Inputs for Mason MDI

# Number of iterations for MDI to perform
num_iter=5000;

# Thinning variable for MDI
thin=1;

# Directory to save output to
save_dir="/home/sdc56/rds/hpc-work/MDI_output/Simulations/Single_dataset/";

# Random seed to use
seed=1;

# Directory to read data from
curr_dir="/home/sdc56/rds/hpc-work/Input_data/Simulations/Single_dataset/Datasets/";

# Scenario to run
scn="base_case";

# Number of clusters allowed in MDI
n_clust=50;

# Simulation number
sim_num=$SLURM_ARRAY_TASK_ID;

# Inference type; one of Consensus (default), Bayesian or Frequentist
model="Consensus";

# For the "base_case" scenario (200 x 20 dataset), 5,000 iterations takes
# about 8 seonds - let's round up.
iter_5k=10;

# Did this get passed to the cluster? 1=True, 0=False
hpc_submit=1;

# Receive parameters from command line
while getopts n:t:o:s:i:k:m:l:h: option; do
    case "${option}" in
        n) num_iter=${OPTARG};;
        t) thin=${OPTARG};;
        o) scn=${OPTARG};;
        s) seed=${OPTARG};;
        k) n_clust=${OPTARG};;
        m) model=${OPTARG};;
        l) iter_5k=${OPTARG};;
        h) hpc_submit=${OPTARG};;
    esac
done

# Output and input homes
input_data="${curr_dir}${scn}/dataset_${sim_num}.csv";
output_dir="${save_dir}${scn}/simulation_${sim_num}/${model}/";

# This is one way to handle different seeds and whatnot
# Instead I have changed the output csv name.
## If Consensus model save to thinning specific directory
#if [ ${model} == "Consensus" ]; then
#  output_dir="${output_dir}/Thin_${thin}/";
#fi


# The output file name
output="${output_dir}${model}N${num_iter}T${thin}Seed${seed}K${n_clust}.csv";

#if [[ -f "$output" ]]; then
#    #echo "$output exist, skipping run."
#    #exit 0
#fi

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
let num_sec=${num_sec}-60*${num_min}-3600*${num_hours};
let num_min=${num_min}-60*${num_hours};

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
###SBATCH --time=${time_str}


#! ############################################################
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:
module load cuda/8.0

#! Full path to application executable: 
application="/home/sdc56/hpc_stuff/scripts/mdipp-1.0.1/mdipp"

options="N ${input_data} -n ${num_iter} -t ${thin} -s ${seed} -c ${n_clust} > ${output}";

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.


cd $workdir

## If not submitting as a job to the cluster, don't do all these things

if [ ${hpc_submit} -eq 1 ]
then

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
#np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.


#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
# CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"

#! (OMP_NUM_THREADS threads will be created):
CMD="{ time $application $options; } 2> ${time_file} ";


###############################################################
### You should not have to change anything below this line ####
###############################################################

# cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
echo "Time allocated: ${SBATCH_TIMELIMIT}";

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"
fi

#! (OMP_NUM_THREADS threads will be created):
# CMD="{ time $application $options; } 2> ${time_file} ";
CMD="$application $options";

eval $CMD 

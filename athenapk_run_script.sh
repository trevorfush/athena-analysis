#!/bin/bash
# ######### SBATCH Lines for Resource Request ##########
 
#SBATCH --time=09:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=2                   # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=128                # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
###SBATCH --gres=gpu:v100:1
#SBATCH --constraint="amd20"
###SBATCH -w amr-089

#SBATCH --mem-per-cpu=8G                     # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name turb_sub_athenapk_rk3_hlld_wenoz     # you can give your job a name for easier identification (same as -J)

#SBATCH -o "runinfo_turb_sub_athenapk_rk3_hlld_wenoz.out"

#SBATCH --mail-type=ALL
#SBATCH --mail-user=fushstep@msu.edu

module purge;

export LANG=en_US.UTF-8

###############################################
#   Setting up simulation numerical methods   #
###############################################

INTEGRATOR=rk3  # rk3, rk2, vl2, rk1
RIEMANN=hlld    # HLLD, HLLE, Roe
RECONST=wenoz     # 2=PLM, 3=PPM, 4=PPM+T

SUBORSUP=sub   # Sub or super (sup) sonic flow

RUNTYPE=cpu    # Hardware being used [gpu or cpu]
NUMCPUS=16     # Number of cpus used for flow analysis

RECOMPILE=1    # 1 = yes, 0 = no... only need to recompile if changing Riemann!

###############################################

SAVEDIR=${HOME}/athenapk/turb-sims/turbulence_10T_${SUBORSUP}_athenapk_${INTEGRATOR}_${RIEMANN}_${RECONST}

# CPU
echo "Running on CPU..."
module load gcccuda/2020a OpenMPI/4.0.3 git HDF5 FFTW CMake

if [ ! -d $SAVEDIR ]
then
    mkdir $SAVEDIR
fi

cd $SAVEDIR

if [ $SUBORSUP = "sub" ]
then 
    # Run with p0 = 1.0 for a Ms ~ 0.5
    echo "Running simulation with Ms = 0.5..."
    srun -n 128 ${HOME}/athenapk/build-cpu/bin/athenaPK -i ${HOME}/athenapk/inputs/turbulence.in parthenon/time/integrator=$INTEGRATOR hydro/reconstruction=$RECONST hydro/riemann=$RIEMANN problem/turbulence/p0=1.0 parthenon/meshblock/nx1=32 parthenon/meshblock/nx2=64 parthenon/meshblock/nx3=64 parthenon/mesh/nghost=4
else
    # Run with p0 = 0.02777 for a Ms ~ 3.0
    echo "Running simulation with Ms = 3.0..."
    srun -n 128 ${HOME}/athenapk/build-cpu/bin/athenaPK -i ${HOME}/athenapk/inputs/turbulence.in parthenon/time/integrator=$INTEGRATOR hydro/reconstruction=$RECONST hydro/riemann=$RIEMANN problem/turbulence/p0=0.02777 parthenon/meshblock/nx1=32 parthenon/meshblock/nx2=64 parthenon/meshblock/nx3=64 parthenon/mesh/nghost=4
fi

    # for num in $(seq -f %05g 0 101); 
    # do 
    #     echo "Running flow analysis on simulation outputs..."
    #     mpirun -n $NUMCPUS python /mnt/home/fushstep/energy-transfer-analysis/run_analysis.py --res 256 --data_path ${num}.athdf --data_type AthenaPP --binning test --type flow --outfile flow_analysis_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}_${num}.hdf5 --eos isothermal --gamma 1.0001  -forced -b;
    # done
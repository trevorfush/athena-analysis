#!/bin/bash
# ######### SBATCH Lines for Resource Request ##########
 
#SBATCH --time=11:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                   # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=128                # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
###SBATCH --gres=gpu:v100:1
#SBATCH --constraint="amd20"
###SBATCH -w amr-089

#SBATCH --mem-per-cpu=6G                     # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name turb_sup_rk4_hlld_ppm     # you can give your job a name for easier identification (same as -J)

#SBATCH -o "runinfo_turb_sup_rk4_hlld_ppm.out"

#SBATCH --mail-type=ALL
#SBATCH --mail-user=fushstep@msu.edu

module purge;

export LANG=en_US.UTF-8

###############################################
#   Setting up simulation numerical methods   #
###############################################

INTEGRATOR=rk4 # rk4, rk3, rk2, vl2
RIEMANN=hlld    # HLLD, HLLE, Roe
RECONST=3      # 2=PLM, 3=PPM, 4=PPM+T

SUBORSUP=sup   # Sub or super (sup) sonic flow

RUNTYPE=cpu    # Hardware being used [gpu or cpu]
NUMCPUS=16     # Number of cpus used for flow analysis

RECOMPILE=1    # 1 = yes, 0 = no... only need to recompile if changing Riemann!

SAVEDIR=/mnt/gs21/scratch/fushstep

###############################################

# GPU
if [ $RUNTYPE = "gpu" ]
then
    echo "Running on GPU..."
    module load gcccuda/2020a OpenMPI/4.0.3 git HDF5 Qt5 FFTW CMake
    export OMPI_CXX=$(pwd)/./kokkos/bin/nvcc_wrapper

    make clean

    ./configure.py --kokkos_arch="Volta70,SKX" --kokkos_devices="Cuda,OpenMP"  --prob=fmturb  --flux=$RIEMANN -hdf5 --hdf5_path=$EBROOTHDF5 -mpi  -b  --nghost=4

    make -j 16

    if [ ! -d ${SAVEDIR}/turb-sims/turbulence_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST} ]
    then
        mkdir ${SAVEDIR}/turb-sims/turbulence_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}
    fi

    cd ${SAVEDIR}/turb-sims/turbulence_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}


    if [ $SUBORSUP = "sub" ]
    then 
        # Run with p0 = 1.0 for a Ms ~ 0.5
        echo "Running simulation with Ms = 0.5..."
        srun -n 1 ${HOME}/kathena/bin/athena -i ${HOME}/kathena/inputs/mhd/athinput.fmturb time/integrator=$INTEGRATOR time/xorder=$RECONST problem/p0=1.0 meshblock/nx1=256 meshblock/nx2=256 meshblock/nx3=256
    else
        # Run with p0 = 0.02777 for a Ms ~ 3.0
        echo "Running simulation with Ms = 3.0..."
        srun -n 1 ${HOME}/kathena/bin/athena -i ${HOME}/kathena/inputs/mhd/athinput.fmturb time/integrator=$INTEGRATOR time/xorder=$RECONST problem/p0=0.05 meshblock/nx1=256 meshblock/nx2=256 meshblock/nx3=256
    fi

    # for num in $(seq -f %05g 0 101); 
    # do 
    #     echo "Running flow analysis on simulation outputs..."
    #     mpirun -n 16 python /mnt/home/fushstep/energy-transfer-analysis/run_analysis.py --res 256 --data_path ${num}.athdf --data_type AthenaPP --binning test --type flow --outfile flow_analysis_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}_${num}.hdf5 --eos isothermal --gamma 1.0001  -forced -b;
    # done

else
    # CPU
    echo "Running on CPU..."
    module load intel/2020a iccifort/2020.4.304 impi  HDF5

    # export OMP_PROC_BIND=true

    if [ $RECOMPILE -eq 1 ]
    then 
        make clean

        ./configure.py --kokkos_arch="AMD" --kokkos_devices="OpenMP"  --prob=fmturb  --flux=$RIEMANN -hdf5 --hdf5_path=$EBROOTHDF5 -mpi  -b --nghost=4 

        make -j 16
    fi

    if [ ! -d ${SAVEDIR}/turb-sims/turbulence_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST} ]
    then
        mkdir ${SAVEDIR}/turb-sims/turbulence_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}
    fi

    cd ${SAVEDIR}/turb-sims/turbulence_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}

    if [ $SUBORSUP = "sub" ]
    then 
        # Run with p0 = 1.0 for a Ms ~ 0.5
        echo "Running simulation with Ms = 0.5..."
        srun -n 128 ${HOME}/kathena/bin/athena -i ${HOME}/kathena/inputs/mhd/athinput.fmturb time/integrator=$INTEGRATOR time/xorder=$RECONST problem/p0=1.0 meshblock/nx1=32 meshblock/nx2=64 meshblock/nx3=64
    else
        # Run with p0 = 0.02777 for a Ms ~ 3.0
        echo "Running simulation with Ms = 3.0..."
        srun -n 128 ${HOME}/kathena/bin/athena -i ${HOME}/kathena/inputs/mhd/athinput.fmturb time/integrator=$INTEGRATOR time/xorder=$RECONST problem/rho0=10.0 problem/p0=0.5 problem/b0=0.1 meshblock/nx1=32 meshblock/nx2=64 meshblock/nx3=64
    fi

    # for num in $(seq -f %05g 0 101); 
    # do 
    #     echo "Running flow analysis on simulation outputs..."
    #     mpirun -n $NUMCPUS python /mnt/home/fushstep/energy-transfer-analysis/run_analysis.py --res 256 --data_path ${num}.athdf --data_type AthenaPP --binning test --type flow --outfile flow_analysis_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}_${num}.hdf5 --eos isothermal --gamma 1.0001  -forced -b;
    # done
fi
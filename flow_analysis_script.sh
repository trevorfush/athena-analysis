#!/bin/bash
# ######### SBATCH Lines for Resource Request ##########
 
#SBATCH --time=08:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                   # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=16                # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --constraint="amd20"
###SBATCH -w amr-089

#SBATCH --mem-per-cpu=6G                     # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name flow_an     # you can give your job a name for easier identification (same as -J)

#SBATCH -o "flow_an.out"

#SBATCH --mail-type=ALL
#SBATCH --mail-user=fushstep@msu.edu

module purge;

export LANG=en_US.UTF-8

###############################################
#   Setting up simulation numerical methods   #
###############################################

# INTEGRATOR=rk4 # rk4, rk3, rk2, vl2
# RIEMANN=roe    # HLLD, HLLE, Roe
# RECONST=4      # 2=PLM, 3=PPM, 4=PPM+T

# SUBORSUP=sub   # Sub or super (sup) sonic flow

# RUNTYPE=cpu    # Hardware being used [gpu or cpu]
# NUMCPUS=16     # Number of cpus used for flow analysis

# RECOMPILE=1    # 1 = yes, 0 = no... only need to recompile if changing Riemann!

###############################################

# CPU
echo "Running on CPU..."
# module load intel/2020a iccifort/2020.4.304 impi  HDF5

SUBORSUP=sub

for RIEMANN in hlle hlld roe;
do
    for RECONST in 2 3 4;
    do
        for INTEGRATOR in vl2 rk2 rk3 rk4;
        do

            if [ -d ${SCRATCH}/turb-sims/turbulence_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST} ]
            then
                cd ${SCRATCH}/turb-sims/turbulence_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}

                echo "Looking at ${SUBORSUP}, ${INTEGRATOR}, ${RIEMANN}, ${RECONST}"

                for num in $(seq -f %05g 0 101); 
                do 

                    if [ ! -f flow_analysis_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}_${num}.hdf5 ]
                    then 
                        echo "Running flow analysis on simulation output ${num}..."
                        mpirun -n 16 python /mnt/home/fushstep/energy-transfer-analysis/run_analysis.py --res 256 --data_path ${num}.athdf --data_type AthenaPP --binning test --type flow --outfile flow_analysis_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}_${num}.hdf5 --eos isothermal --gamma 1.0001  -forced -b;
                    fi
                    
                done

            fi

        done
    done
done

# cd ${SCRATCH}/turb-sims/turbulence_10T_${SUBORSUP}_${INTEGRATOR}_${RIEMANN}_${RECONST}

# if [ $SUBORSUP = "sub" ]
# then 
#     # Run with p0 = 1.0 for a Ms ~ 0.5
#     echo "Running simulation with Ms = 0.5..."
#     srun -n 128 ${HOME}/kathena/bin/athena -i ${HOME}/kathena/inputs/mhd/athinput.fmturb time/integrator=$INTEGRATOR time/xorder=$RECONST problem/p0=1.0 meshblock/nx1=32 meshblock/nx2=64 meshblock/nx3=64
# else
#     # Run with p0 = 0.02777 for a Ms ~ 3.0
#     echo "Running simulation with Ms = 3.0..."
#     srun -n 128 ${HOME}/kathena/bin/athena -i ${HOME}/kathena/inputs/mhd/athinput.fmturb time/integrator=$INTEGRATOR time/xorder=$RECONST problem/p0=0.02777 meshblock/nx1=32 meshblock/nx2=64 meshblock/nx3=64
# fi


# athena-analysis
### Analysis pipeline for turbulence simulations for calculating norms and divergence

## To run

First, all the simulation directories need to be located in THE SAME DIRECTORY. For example, having a directory `turb-sims` whose contents are ```sim1 sim2 sim3``` each containing either K-Athena turbulence outputs or Athena-PK turbulence outputs. Within this `turb-sims` directory, for the simulation with the combination of numerical methods being used as a reference solution, the directory should be renamed with `REF_` preceeding the directory name. For example, with the three simulations in the above sample directory if `sim1` is the reference simulation, then it should be renamed `REF_sim1`.

There needs to be a reference solution for the L1, L2, and Linf norms to be calculated.

Run python analysismain.py with the arguments...

| Argument | Function |
| ---------|--------- |
| `-p`, `--prob` | Select the problem to be solved: (sedov_noB, sedov_Bx, sedov_By, brio_wu, kh_noB, current_sheet, orszag_tang, turbulence) |
| `--simloc` | Simulation directory where all simulation runs are located |
| `--outfile` | Output data name for analysis quantities |
| `--saveimages` | Whether or not to save the plots associated withpipeline. If turned on, plots will be saved, otherwise just hdf5 data is saved. `turbulence` runs with `--saveimages` on has not been tested. |

## Directory contents

In addition to the norm and divergence calculations, scripts are included that ran the simulation for K-Athena and Athena-PK, as well as scripts to run the flow analysis on each of the turbulence simulations in a directory. To run Athena-PK (K-Athena) simulations, the `athenapk_run_script.sh` (`athenapk_run_script.sh`) can be run after changing the location of the `athenapk` (`kathena`) directory. To run flow analysis, browse through the `flow_analysis_script.sh` script and adjust needed parameters and run.

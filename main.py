import os
import sys
import numpy as np 
import glob
import time
import subprocess
import yaml
# sys.path.append("..")

from src.sim import Runsim
from src.analysis import Analysis

import argparse

#######################################

# RUNNAME IS THE NAME IN THE CONFIG FILE!!!!
# PROBLEM NAME IS THE NAME OF THE CORRESPONDING INPUT FILE IN ATHENA-PK
# PROBLEM ID IS THE PROBLEM ID FROM THE INPUT FILE FOR ATHENA PK

sedov_noB_sim_params = {
                    "runname"        : "sedov_noB",
                    "problem_name"   : "sedov.in",
                    "problem_id"     : "sedov",
                    "riemann"        : None,
                    "reconstruction" : None,
                    "reference_sol"  : None
                   }

brio_wu_sim_params = {
                   "runname"        : "brio_wu",
                   "problem_name"   : "sod.in",
                   "problem_id"     : "sod",
                   "riemann"        : None,
                   "reconstruction" : None,
                   "reference_sol"  : None
                  }

sedov_Bx_sim_params = {
                    "runname"        : "sedov_Bx",
                    "problem_name"   : "sedov.in",
                    "problem_id"     : "sedov",
                    "riemann"        : None,
                    "reconstruction" : None,
                    "reference_sol"  : None
                  }

sedov_By_sim_params = {
                    "runname"        : "sedov_By",
                    "problem_name"   : "sedov.in",
                    "problem_id"     : "sedov",
                    "riemann"        : None,
                    "reconstruction" : None,
                    "reference_sol"  : None
                  }

kh_noB_sim_params = {
                    "runname"        : "kh_noB",
                    "problem_name"   : "kh-shear-lecoanet_2d.in",
                    "problem_id"     : "kh",
                    "riemann"        : None,
                    "reconstruction" : None,
                    "reference_sol"  : None
                }

current_sheet_sim_params = {
                        "runname"        : "current_sheet",
                        "problem_name"   : "current_sheet.in",
                        "problem_id"     : "current_sheet",
                        "riemann"        : None,
                        "reconstruction" : None,
                        "reference_sol"  : None
                       }

orszag_tang_sim_params = {
                        "runname"        : "orszag_tang",
                        "problem_name"   : "orszag_tang.in",
                        "problem_id"     : "orszag_tang",
                        "riemann"        : None,
                        "reconstruction" : None,
                        "reference_sol"  : None
                     }

#######################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-p","--prob", 
                        help="Select the problem to be solved: (sedov_noB, sedov_Bx, sedov_By, brio_wu, kh_noB, current_sheet, orszag_tang, turbulence)")
    
    parser.add_argument("-c","--code", 
                        help="Select the code used to run sim: (athenapk, kathena)")

    parser.add_argument(
        "--riemann", help="Riemann solver to run if not running all possible combinations"
    )

    parser.add_argument(
        "--reconstruction", help="Reconstruction method if not running all possible combinations"
    )

    parser.add_argument(
        "-a", "--all", help="Runs given simulation with all possible combinations of numerical methods for the given code",
        action="store_true"
    )

    parser.add_argument(
        "--reference", help="Determines whether the solution is a reference (flag on) or not (flag off)", action="store_true"
    )

    args = parser.parse_args()

    print(f"args.prob = {args.prob}")
    print(f"args.code = {args.reference}")

    runname = args.prob

    if runname == "sedov_noB":
        params = sedov_noB_sim_params
    elif runname == "sedov_Bx":
        params = sedov_Bx_sim_params
    elif runname == "sedov_By":
        params = sedov_By_sim_params
    elif runname == "brio_wu":
        params = brio_wu_sim_params
    elif runname == "kh_noB":
        params = kh_noB_sim_params
    elif runname == "current_sheet":
        params = current_sheet_sim_params
    elif runname == "orszag_tang":
        params = orszag_tang_sim_params
    else:
        print("[ERROR] Please use a valid implemented runname!!")

    ## Overriding params to run on different combinations of num. meth.
    if args.all == True:

        assert args.reference == False, "-a flag only loops through non-reference solutions, reference solutions should be generated separately"

        for riemann in ["hlld", "hlle"]:
            for recon in ["plm","ppm"]:

                params["riemann"] = riemann
                params["reconstruction"] = recon

                test = Runsim(args, params)
        
                print("Runsim object created, running sim!")
                test.runsim()
                print("Finished running sim!")
    
    else:

        test = Runsim(args, params)

        print("Runsim object created, running sim!")
        test.runsim()
        print("Finished running sim!")
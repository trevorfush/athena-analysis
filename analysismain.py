import os
import string
import sys
import numpy as np 
import glob
import time
import subprocess
import yaml
sys.path.append("..")

from src.sim import Runsim
from src.parallelanalysis import Analysis

import argparse

#######################################

# RUNNAME IS THE NAME IN THE CONFIG FILE!!!!
# PROBLEM NAME IS THE NAME OF THE CORRESPONDING INPUT FILE IN ATHENA-PK
# PROBLEM ID IS THE PROBLEM ID FROM THE INPUT FILE FOR ATHENA PK

sedov_noB_params = {
                    "runname"        : "sedov_noB",
                    "problem_name"   : "sedov.in",
                    "problem_id"     : "sedov",
                    "reconstruction" : "wenoz",
                    "riemann"        : "hlle",
                    "reference_sol"  : False
                   }

brio_wu_params = {
                   "runname"        : "brio_wu",
                   "problem_name"   : "sod.in",
                   "problem_id"     : "sod",
                   "reconstruction" : "plm",
                   "riemann"        : "hlle",
                   "reference_sol"  : False
                  }

sedov_Bx_params = {
                    "runname"        : "sedov_Bx",
                    "problem_name"   : "sedov.in",
                    "problem_id"     : "sedov",
                    "reconstruction" : "plm",
                    "riemann"        : "hlle",
                    "reference_sol"  : False
                  }

sedov_By_params = {
                    "runname"        : "sedov_By",
                    "problem_name"   : "sedov.in",
                    "problem_id"     : "sedov",
                    "reconstruction" : "plm",
                    "riemann"        : "hlle",
                    "reference_sol"  : False
                  }

kh_noB_params = {
                    "runname"        : "kh_noB",
                    "problem_name"   : "kh-shear-lecoanet_2d.in",
                    "problem_id"     : "kh",
                    "reconstruction" : "plm",
                    "riemann"        : "hlle",
                    "reference_sol"  : False
                }

current_sheet_params = {
                        "runname"        : "current_sheet",
                        "problem_name"   : "current_sheet.in",
                        "problem_id"     : "current_sheet",
                        "reconstruction" : "plm",
                        "riemann"        : "hlle",
                        "reference_sol"  : False
                       }

orszag_tang_params = {
                        "runname"        : "orszag_tang",
                        "problem_name"   : "orszag_tang.in",
                        "problem_id"     : "orszag_tang",
                        "reconstruction" : "plm",
                        "riemann"        : "hlle",
                        "reference_sol"  : False
                     }

turbulence_params  = {
    "runname" : "turbulence",
    "problem_name" : "athinput.fmturb",
    "problem_id"   : "Turb",
    "reconstruction" : "ppm",
    "riemann"      : "hlle",
    "reference_sol" : False
}

#######################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-p","--prob", 
                        help="Select the problem to be solved: (sedov_noB, sedov_Bx, sedov_By, brio_wu, kh_noB, current_sheet, orszag_tang, turbulence)")

    parser.add_argument(
        "--simloc",
        help="Simulation directory where all simulation runs are located",
        type=str
    )

    parser.add_argument(
        "--outfile",
        help="Output data name for analysis quantities",
        type=str 
    )

    parser.add_argument(
        "--saveimages",
        help="Whether or not to save the plots associated with pipeline. If turned on, plots will be saved, otherwise just hdf5 data is saved",
        action="store_true"
    )

    args = parser.parse_args()

    print(f"args.prob = {args.prob}")
    # print(f"args.code = {args.code}")

    runname = args.prob

    if runname == "sedov_noB":
        params = sedov_noB_params
    elif runname == "sedov_Bx":
        params = sedov_Bx_params
    elif runname == "sedov_By":
        params = sedov_By_params
    elif runname == "brio_wu":
        params = brio_wu_params
    elif runname == "kh_noB":
        params = kh_noB_params
    elif runname == "current_sheet":
        params = current_sheet_params
    elif runname == "orszag_tang":
        params = orszag_tang_params
    elif runname == "turbulence":
        params = turbulence_params
    else:
        print("[ERROR] Please use a valid implemented runname!!")

    test = Analysis(params, args.simloc, args.outfile, args.saveimages)

    test.Analyze()
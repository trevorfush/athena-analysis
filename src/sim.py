import os
import sys
import numpy as np 
import glob
import time
import subprocess
import yaml
sys.path.append(os.path.join("..", ".."))


class Runsim():

    def __init__(self, params):
        
        # SHOULD MAKE FILES WITH PARAMETERS FOR EACH OF THE SIMULATIONS
        # TO BE RUN, THEN READ THOSE TO INCLUDE IN THIS CLASS
        
        # RUNNAME IS THE NAME IN THE CONFIG FILE!!!!
        self.runname   = params["runname"]

        # PROBLEM NAME IS THE NAME OF THE CORRESPONDING INPUT FILE IN ATHENA-PK
        self.inputfile = params["problem_name"]

        # PROBLEM ID IS THE PROBLEM ID FROM THE INPUT FILE FOR ATHENA PK
        self.problem_id = params["problem_id"]

        self.recon     = params["reconstruction"]
        self.riemann   = params["riemann"]

        self.reference = params["reference_sol"]

        self.simparms  = {}

        self.path      = "data"

        self.athenapk  = "../../athenapk"

        if self.reference == True:
            self.savename = f"REF_{self.runname}"
        else:
            self.savename = self.runname

    def load_params(self):

        with open("simconfig.yml", "r") as f:
            config = yaml.safe_load(f)
            f.close()

        self.simparams = config[self.runname]

    def checkprogress(self):
        """
        Checks progress of current simulation of Runsim object. Returns
        a 1x3 array:

        current_progress[0] --> Has simulation started running
        current_progress[1] --> If started running, is it completed
        current_progress[2] --> If not complete, where did it leave off 
        """

        current_progress = [False, False, None]

        direxist = os.path.isdir(f"{self.path}/{self.savename}")
        fileexist = os.path.isfile(f"{self.path}/{self.savename}/FINISHED.txt")

        print(direxist, fileexist)

        if (direxist == True and len(os.listdir(f"{self.path}/{self.savename}")) > 0):
            current_progress[0] = True
        if direxist == False:
            current_progress[0] = False

        if fileexist == True:
            current_progress[1] = True
            current_progress[2] = "DONE"
        if  (fileexist == False and direxist == True):
            current_progress[1] = False
            current_progress[2] = max(glob.glob(f"{self.path}/{self.savename}/*.?????.phdf"))
        if (fileexist == False and direxist == False):
            current_progress[1] = False
            current_progress[2] = None

        return current_progress

    def genFINISHED(self, runtime):

        last_ran = max(glob.glob("./*.?????.phdf"))
        
        # Calculating the last POSSIBLE output to exist based on sim params
        tlim = self.simparams["tlim"]
        dt   = self.simparams["dt"]

        lastnum = str(int(tlim/dt)).zfill(5)

        maxnum_ran = 0
        for (i, item) in enumerate(last_ran.split(".")):
            for subitem  in item.split():
                if (subitem.isdigit()):
                    maxnum_ran = last_ran.split(".")[i]

        if maxnum_ran == lastnum:
        
            with open(f"FINISHED.txt", "w+") as f:

                f.write(f"Finished run with following params: \nRUNTIME = {runtime}\nRECONST = {self.recon}\nRIEMANN = {self.riemann}\nREFERENCE = {self.reference}\n")
                f.close()

        else:
            print("[ERROR] Did not finish running simulation")
            pass

    def runsim(self):
        """
        Runs current simulation based on current progress.
        """

        # Check if simulation has already ran
        cp = self.checkprogress()
        
        # Load config file parameters
        self.load_params()

        #######################################################################
        #  Set up the command line arguments to initialize problem uniformly  #
        #######################################################################
        param_str = ""
        keys = list(self.simparams.keys())
        vals = list(self.simparams.values())

        for i in range(len(self.simparams)):
            if (keys[i] != "fluid") and (keys[i] != "tlim") and (keys[i] != "dt"):
                param_str += f"problem/{self.problem_id}/{keys[i]}={vals[i]}"
            elif (keys[i] != "tlim") and (keys[i] != "dt"):
                param_str += f"hydro/fluid={vals[i]}"
            else:
                pass
        #######################################################################

        # RUN IF The SIMULATION HASN'T BEEN RUN YET
        if cp[0] == False:
            
            # print(here)
            os.mkdir(f"{self.path}/{self.savename}_{self.riemann}_{self.recon}")
            os.chdir(f"{self.path}/{self.savename}_{self.riemann}_{self.recon}")

            # print(self.riemann, self.recon)
            start = time.time()
            os.system(f"../{self.athenapk}/build-host/bin/athenaPK -i ../{self.athenapk}/inputs/{self.inputfile} {param_str} hydro/reconstruction={self.recon} hydro/riemann={self.riemann} parthenon/mesh/nghost=3")
            end = time.time()
            
            rt = end-start

            ## NEED TO CHECK IF IT ACTUALLY FINISHED HERE
            self.genFINISHED(runtime=rt)
        
        # RUN PICKING UP WHERE THE SIMULATION LEFT OFF
        elif (cp[0] == True and cp[1] == False):

            last = cp[2]

            # NEED TO FILL THIS IN
            os.chdir(f"{self.path}/{self.savename}")
            os.system(f"../{self.athenapk}/build-host/bin/athenaPK -i ../{self.athenapk}/inputs/{self.inputfile} -r {last} hydro/reconstruction={self.recon} hydro/riemann={self.riemann} parthenon/mesh/nghost=3")

        else:

            print("Something went weird, might want to double check something!")


# if __name__ == "__main__":

    # RUNNAME IS THE NAME IN THE CONFIG FILE!!!!
    # PROBLEM NAME IS THE NAME OF THE CORRESPONDING INPUT FILE IN ATHENA-PK
    # PROBLEM ID IS THE PROBLEM ID FROM THE INPUT FILE FOR ATHENA PK

    # sedov_noB_params = {
    #                     "runname"        : "sedov_noB",
    #                     "problem_name"   : "sedov.in",
    #                     "problem_id"     : "sedov",
    #                     "reconstruction" : "plm",
    #                     "riemann"        : "hlle",
    #                     "reference_sol"  : False
    #                    }

    # brio_wu_params = {
    #                   "runname"        : "brio_wu",
    #                   "problem_name"   : "sod.in",
    #                   "problem_id"     : "sod",
    #                   "reconstruction" : "plm",
    #                   "riemann"        : "hlle",
    #                   "reference_sol"  : False
    #                  }

    # sedov_Bx_params = {
    #                    "runname"        : "sedov_Bx",
    #                    "problem_name"   : "sedov.in",
    #                    "problem_id"     : "sedov",
    #                    "reconstruction" : "plm",
    #                    "riemann"        : "hlle",
    #                    "reference_sol"  : False
    #                   }

    # sedov_By_params = {
    #                    "runname"        : "sedov_By",
    #                    "problem_name"   : "sedov.in",
    #                    "problem_id"     : "sedov",
    #                    "reconstruction" : "plm",
    #                    "riemann"        : "hlle",
    #                    "reference_sol"  : False
    #                   }

    # kh_noB_params = {
    #                  "runname"        : "kh_noB",
    #                  "problem_name"   : "kh-shear-lecoanet_2d.in",
    #                  "problem_id"     : "kh",
    #                  "reconstruction" : "plm",
    #                  "riemann"        : "hlle",
    #                  "reference_sol"  : False
    #                 }

    # current_sheet_params = {
    #                         "runname"        : "current_sheet",
    #                         "problem_name"   : "current_sheet.in",
    #                         "problem_id"     : "current_sheet",
    #                         "reconstruction" : "plm",
    #                         "riemann"        : "hlle",
    #                         "reference_sol"  : False
    #                        }

    # orszag_tang_params = {
    #                       "runname"        : "orszag_tang",
    #                       "problem_name"   : "orszag_tang.in",
    #                       "problem_id"     : "orszag_tang",
    #                       "reconstruction" : "plm",
    #                       "riemann"        : "hlle",
    #                       "reference_sol"  : False
    #                      }

    # test = Runsim(params)
    
    # print("Runsim object created, running sim!")
    # test.runsim()
    # print("Finished running sim!")
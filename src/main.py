import os
import sys
import numpy as np 
import glob
import time
import subprocess
sys.path.append(os.path.join("..", ".."))


class Runsim():

    def __init__(self, params):
        
        # SHOULD MAKE FILES WITH PARAMETERS FOR EACH OF THE SIMULATIONS
        # TO BE RUN, THEN READ THOSE TO INCLUDE IN THIS CLASS

        self.runname   = params["runname"]
        self.inputfile = params["problem_name"]

        self.recon     = params["reconstruction"]
        self.riemann   = params["riemann"]

        self.reference = params["reference_sol"]

        self.max_out   = params["last_output"]

        self.path      = "../data"

        self.athenapk  = "../../athenapk"

        if self.reference == True:
            self.savename = f"REF_{self.runname}"
        else:
            self.savename = self.runname


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

        last_ran = max(glob.glob(f"{self.path}/{self.savename}/*.?????.phdf"))

        if self.max_out == last_ran:
        
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

        cp = self.checkprogress()

        # print(type(cp[0]))

        if cp[0] == False:
            
            # print(here)
            os.mkdir(f"{self.path}/{self.savename}")
            os.chdir(f"{self.path}/{self.savename}")

            # print(self.riemann, self.recon)
            start = time.time()
            os.system(f"../{self.athenapk}/build-host/bin/athenaPK -i ../{self.athenapk}/inputs/{self.inputfile} hydro/reconstruction={self.recon} hydro/riemann={self.riemann} parthenon/mesh/nghost=3")
            end = time.time()
            
            rt = end-start

            ## NEED TO CHECK IF IT ACTUALLY FINISHED HERE
            self.genFINISHED(runtime=rt)

        elif (cp[0] == True and cp[1] == False):

            last = cp[2]

            # NEED TO FILL THIS IN
            os.chdir(f"{self.path}/{self.savename}")
            os.system(f"../{self.athenapk}/build-host/bin/athenaPK -i ../{self.athenapk}/inputs/{self.inputfile} -r {last} hydro/reconstruction={self.recon} hydro/riemann={self.riemann} parthenon/mesh/nghost=3")

        else:

            print("Something went weird, might want to double check something!")


if __name__ == "__main__":

    params = {"runname"        : "sedov_restart_test",
              "problem_name"   : "sedov.in",
              "reconstruction" : "plm",
              "riemann"        : "hlle",
              "reference_sol"  : False,
              "last_output"    : "parthenon.prim.00020.phdf"}

    test = Runsim(params)
    
    print("Runsim object created, running sim!")
    test.runsim()
    print("Finished running sim!")
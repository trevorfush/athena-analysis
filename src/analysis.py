import numpy as np
import yt
import matplotlib.pyplot as plt 
import matplotlib as mpl 
mpl.use("Agg")
import os
import sys
import glob
import time


class Analysis():

    def __init__(self, params):

        self.runname   = params["runname"]
        self.inputfile = params["problem_name"]

        self.recon     = params["reconstruction"]
        self.riemann   = params["riemann"]

        self.reference = params["reference_sol"]

        self.max_out   = params["last_output"]

        self.path      = "../results"
        self.simpath   = "../data"

        self.athenapk  = "../../athenapk"

        if self.reference == True:
            self.savename = f"REF_{self.runname}"
        else:
            self.savename = self.runname

        self.plotsavepreamble = f"{self.runname.split(".")[0]}_{self.riemann}_{self.recon}"

        self.dslist = glob.glob(f"{self.simpath}/{self.savename}/*.?????.phdf"
        self.maxds  = yt.load(max(self.dslist))


    def checkSimProgress(self):
        """
        Checks progress of current simulation of Runsim object. Returns
        a 1x3 array:

        current_progress[0] --> Has simulation started running
        current_progress[1] --> If started running, is it completed
        current_progress[2] --> If not complete, where did it leave off 
        """

        current_progress = [False, False, None]

        direxist = os.path.isdir(f"{self.simpath}/{self.savename}")
        fileexist = os.path.isfile(f"{self.simpath}/{self.savename}/FINISHED.txt")

        print(direxist, fileexist)

        if (direxist == True and len(os.listdir(f"{self.simpath}/{self.savename}")) > 0):
            current_progress[0] = True
        if direxist == False:
            current_progress[0] = False

        if fileexist == True:
            current_progress[1] = True
            current_progress[2] = "DONE"
        if  (fileexist == False and direxist == True):
            current_progress[1] = False
            current_progress[2] = max(glob.glob(f"{self.simpath}/{self.savename}/*.?????.phdf"))
        if (fileexist == False and direxist == False):
            current_progress[1] = False
            current_progress[2] = None

        return current_progress

    def getFinalSlices(self):

        start = time.time()

        for field in self.maxds.field_list:
            
            title = f"{self.plotsavepreamble}_{field}"

            p = yt.SlicePlot(self.maxds, fields=field, axis="z")
            p.annotate_timestamp(time_format="t = {time:.3f}", text_args={'color':'white',
                                'horizontalalignment':'center', 'verticalalignment':'top', 
                                "size":20})

            p.annotate_title(title)
            p.save(name=f"{title}.png")

        end = time.time()

        with open(f"SLICES_FINISHED.txt", "w+") as f:

                f.write(f"Finished with slices in {end - start} sec.")
                f.close()

    def checkFinalSlices(self):

        return os.path.isfile(f"{self.path}/{self.savename}/SLICES_FINISHED.txt")

    def getNorms(self, ):

    def checkNorms(self):

        return os.path.isfile(f"{self.path}/{self.savename}/NORMS_FINISHED.txt")

    def getFinalErrors(self, ):

    def checkFinalErrors(self):

        return os.path.isfile(f"{self.path}/{self.savename}/ERRORS_FINISHED.txt")

    def checkForReference(self):

        sims  = ["sod", "sedov", "kh", "current_sheet", "turbulence", "orszag_tang"]
        exist = [False for i in range(len(sims))]

        for (i, simtype) in enumerate(sims):

            if os.path.isfile(f"{self.simpath}/REF_{simtype}/FINISHED.txt"):
                exist[i] = True

        if False not in exist:
            return {"success":0, "exist":exist}
        else:
            return {"success":1, "exist":exist}


    def getWebpage(self):
        """
        Generates HTML file with all of the plots
        """

    def Analyze(self):
        """
        Runs analysis on given simulation and tracks progress.
        """

        # Double check that simulation is finished

        cp = self.checkprogress()
        assert cp[1] == True "[ERROR] Simulation not finished running!!!"

        # Create a directory for the given simulation in the results directory

        os.mkdir(f"{self.path}/{self.savename}")
        os.chdir(f"{self.path}/{self.savename}")

        # Check if slices complete. If not, run getFinalSlices

        slice_comp = self.checkFinalSlices()
        if slice_comp == True:
            pass
        else:
            self.getFinalSlices()

        # Check if norms plot complete. If not, run getNorms



        # Check if final errors compelte. If not, run getFinalErrors



        # Generate HTML file/webpage



        return 0




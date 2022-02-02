import numpy as np
import yt
import matplotlib.pyplot as plt 
import matplotlib as mpl 
mpl.use("Agg")
import os
import sys
import glob
import time

yt.set_log_level(50)

class Analysis():

    def __init__(self, params):


        # RUNNAME IS THE NAME IN THE CONFIG FILE!!!!
        self.runname   = params["runname"]

        # PROBLEM NAME IS THE NAME OF THE CORRESPONDING INPUT FILE IN ATHENA-PK
        self.inputfile = params["problem_name"]

        # PROBLEM ID IS THE PROBLEM ID FROM THE INPUT FILE FOR ATHENA-PK
        self.problem_id = params["problem_id"]


        self.recon     = params["reconstruction"]
        self.riemann   = params["riemann"]

        self.reference = params["reference_sol"]

        self.path      = "results"
        self.simpath   = "data"

        self.athenapk  = "../athenapk"

        if self.reference == True:
            self.savename = f"REF_{self.runname}"
        else:
            self.savename = self.runname

        
        self.plotsavepreamble = f"{self.runname}_{self.riemann}_{self.recon}"

        ###############################################################################

        # This needs to change!!!!!!!!!!!!!!

        # self.dslist = glob.glob(f"{self.simpath}/{self.savename}/*.?????.phdf")
        # self.maxds  = yt.load(max(self.dslist))

        ###############################################################################
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

    def getNormVals(self, output, riemann, recon, field):

        # Getting reference riemann and reconstruction
        reference_strs = glob.glob(f"{self.simpath}/REF_{self.runname}_*")
        reference_str  = reference_strs[0]

        ref_riemann    = reference_str.split("_")[-2]
        ref_recon      = reference_str.split("_")[-1]

        ds     = yt.load(output)
        ref_ds = yt.load(f"{self.simpath}/REF_{self.runname}_{ref_riemann}_{ref_recon}/{output.split('/')[-1]}")
        t = np.float64(ds.current_time.value)

        ad  = ds.all_data()
        arr = ad[field]
        data = arr.to_ndarray().reshape(ds.domain_dimensions)

        ad  = ref_ds.all_data()
        arr = ad[field]
        ref_data = arr.to_ndarray().reshape(ref_ds.domain_dimensions)

        N = 1
        for dim in data.shape:
            N *= dim

        e_ij = np.abs(np.abs(ref_data) - np.abs(data))

        l1 = 1/N * (e_ij.sum())
        l2 = np.sqrt(1/N * ((e_ij**2).sum()))
        linf = np.amax(e_ij)

        return l1, l2, linf, t


    def getFinalSlices(self):

        start = time.time()

        sims = glob.glob(f"{self.simpath}/{self.savename}_*")

        for sim in sims:
            print(f"Making slices for {sim} ...")
            riemann = sim.split("_")[-2]
            reconst = sim.split("_")[-1]

            output = max(sorted(glob.glob(f"{sim}/*.phdf")))

            ds = yt.load(output)

            for field in ds.field_list:
                
                title = f"{self.runname}_{riemann}_{reconst}_{field[-1]}"

                p = yt.SlicePlot(ds, "z", field=field)
                p.annotate_timestamp(time_format="t = {time:.3f}", text_args={'color':'white',
                                    'horizontalalignment':'center', 'verticalalignment':'top', 
                                    "size":20})

                p.annotate_title(title)
                p.save(name=f"{self.path}/{self.runname}/{title}.png")

        end = time.time()

        with open(f"{self.path}/{self.runname}/SLICES_FINISHED.txt", "w+") as f:

                f.write(f"Finished with slices in {end - start} sec. \n")
                f.close()

    def checkFinalSlices(self):

        return os.path.isfile(f"{self.path}/{self.savename}/SLICES_FINISHED.txt")

    def getNorms(self):

        start = time.time()

        sims = glob.glob(f"{self.simpath}/{self.savename}_*")

        # NEED TO CHANGE THIS TO DO ALL OF THE FIELDS (MAYBE???)

        field = "Density"

        riemanns = []
        recons   = []

        l1s_bysim = []
        l2s_bysim = []
        linfs_bysim = []
        ts_bysim    = []

        for sim in sims:
            riemann = sim.split("_")[-2]
            reconst = sim.split("_")[-1]

            riemanns.append(riemann)
            recons.append(reconst)

            outputs = sorted(glob.glob(f"{sim}/*.phdf"))

            l1s = []
            l2s = []
            linfs = []
            ts    = []

            for output in outputs:
                print(output)
                l1, l2, linf, t = self.getNormVals(output, riemann, reconst, field)

                l1s.append(l1)
                l2s.append(l2)
                linfs.append(linf)
                ts.append(t)

            l1s_bysim.append(l1s)
            l2s_bysim.append(l2s)
            linfs_bysim.append(linfs)
            ts_bysim.append(ts)

        fig, ax = plt.subplots(3, 1, sharex=True, figsize=(8,10))

        cols = plt.cm.cool(np.linspace(0,1,len(sims)))

        for i in range(len(l1s_bysim)):

            ax[0].plot(ts_bysim[i], l1s_bysim[i], color=cols[i], label=f"{riemanns[i]}_{recons[i]}")

            ax[1].plot(ts_bysim[i], l2s_bysim[i], color=cols[i], label=f"{riemanns[i]}_{recons[i]}")
        
            ax[2].plot(ts_bysim[i], linfs_bysim[i], color=cols[i], label=f"{riemanns[i]}_{recons[i]}")

        ax[0].set_ylabel(r"$L_1$", rotation=0)

        ax[1].set_ylabel(r"$L_2$", rotation=0)

        ax[2].set_ylabel(r"$L_{\infty}$", rotation=0)
        ax[2].set_xlabel(r"Time [code time]")

        for i in range(3):
            ax[i].legend()
            ax[i].grid()

        savetitle = f"{self.path}/{self.savename}/norms_vs_time.png"

        plt.savefig(savetitle)

        plt.close()

        end = time.time()

        with open(f"{self.path}/{self.runname}/NORMS_FINISHED.txt", "w+") as f:

                f.write(f"Finished with slices in {end - start} sec. \n")
                f.close()



    def checkNorms(self):

        return os.path.isfile(f"{self.path}/{self.savename}/NORMS_FINISHED.txt")

    def getFinalErrors(self, ):

        return

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
        #######################################################################
        
        # THIS NEEDS TO CHANGE TOO, NOT JUST ONE SIMULATION IS BEING ANALYZED

        # print("Checking Progress")
        # cp = self.checkSimProgress()
        # assert cp[1] == True, "[ERROR] Simulation not finished running!!!"
        # print("Finished checking progress, proceeding")
        #######################################################################

        # Create a directory for the given simulation in the results directory

        if os.path.isdir(f"{self.path}/{self.savename}") != True:
            os.mkdir(f"{self.path}/{self.savename}")

        # os.chdir(f"{self.path}/{self.savename}")

        # Check if slices complete. If not, run getFinalSlices

        print("Getting slices")
        slice_comp = self.checkFinalSlices()
        if slice_comp == True:
            print("Already finished getting slices!")
            pass
        else:
            print("RUNNING SLICES")
            self.getFinalSlices()
            print("Finished getting slices")

        # Check if norms plot complete. If not, run getNorms

        print("Getting norms")
        norm_comp = self.checkNorms()
        if norm_comp == True:
            print("Already finished getting norms!")
            pass
        else:
            print("RUNNING NORMS")
            self.getNorms()
            print("Finished getting norms")

        # Check if final errors compelte. If not, run getFinalErrors



        # Generate HTML file/webpage



        return 0




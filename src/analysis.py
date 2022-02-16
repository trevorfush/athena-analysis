import numpy as np
import yt
import matplotlib.pyplot as plt 
import matplotlib as mpl 
mpl.use("Agg")
import os
import sys
import glob
import time
from tqdm import tqdm
import yaml
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
        
        self.simparams = 0

        ###############################################################################

        # This needs to change!!!!!!!!!!!!!!

        # self.dslist = glob.glob(f"{self.simpath}/{self.savename}/*.?????.phdf")
        # self.maxds  = yt.load(max(self.dslist))

        ###############################################################################

    def load_params(self):

        with open("simconfig.yml", "r") as f:
            config = yaml.safe_load(f)
            f.close()

        self.simparams = config[self.runname]

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

    def getNormVals(self, ind, ds, refds, l1_array, l2_array, linf_array, t_array):

        # Getting reference riemann and reconstruction
        # reference_strs = glob.glob(f"{self.simpath}/REF_{self.runname}_*")
        # reference_str  = reference_strs[0]

        # ref_riemann    = reference_str.split("_")[-2]
        # ref_recon      = reference_str.split("_")[-1]

        # ds     = yt.load(output)
        # ref_ds = yt.load(f"{self.simpath}/REF_{self.runname}_{ref_riemann}_{ref_recon}/{output.split('/')[-1]}")
        
        t = np.float64(ds.current_time.value)

        ad  = ds.all_data()
        refad  = refds.all_data()

        for (j, field) in enumerate(ds.field_list):

            arr = ad[field]
            data = arr.to_ndarray().reshape(ds.domain_dimensions)

            refarr = refad[field]
            ref_data = refarr.to_ndarray().reshape(refds.domain_dimensions)

            N = 1
            for dim in data.shape:
                N *= dim

            # e_ij = np.abs(np.abs(ref_data) - np.abs(data)) / np.abs(ref_data)
            # NEED TO GET THIS TO WORK WITH SMALL VALUES SOMEHOW

            e_ij = np.abs(np.abs(ref_data) - np.abs(data))

            l1 = 1/N * (e_ij.sum())
            l2 = np.sqrt(1/N * ((e_ij**2).sum()))
            linf = np.amax(e_ij)

            l1_array[j, ind] = l1
            l2_array[j, ind] = l2
            linf_array[j, ind] = linf
            t_array[ind] = t

    def getDifferencePlots(self, ind, ds, refds, t_array, riemann, reconst, ref_riemann, ref_reconst):

        ad = ds.all_data()
        ref_ad = refds.all_data()

        num_fields = len(ds.field_list)
        fig, ax = plt.subplots(3, num_fields, figsize=(15,9))

        for (j,field) in enumerate(ds.field_list):

            arr = ad[field]
            ref_arr = ref_ad[field]

            data = arr.to_ndarray().reshape(refds.domain_dimensions)
            ref_data = ref_arr.to_ndarray().reshape(refds.domain_dimensions)

            diff = ref_data - data

            im1 = ax[0,j].imshow(ref_data[:,:,0])
            im2 = ax[1,j].imshow(data[:,:,0])
            im3 = ax[2,j].imshow(diff[:,:,0], cmap="jet")

            divider1 = make_axes_locatable(ax[0,j])
            cax1 = divider1.append_axes("right",size="5%",pad=0.05)
            cbar1 = plt.colorbar(im1, cax=cax1)

            divider2 = make_axes_locatable(ax[1,j])
            cax2 = divider2.append_axes("right",size="5%",pad=0.05)
            cbar2 = plt.colorbar(im2, cax=cax2)

            divider3 = make_axes_locatable(ax[2,j])
            cax3 = divider3.append_axes("right",size="5%",pad=0.05)
            cbar3 = plt.colorbar(im3, cax=cax3)

            ax[0,j].set_title(field[-1])

        ax[0,0].set_ylabel("Reference")
        ax[1,0].set_ylabel("Comparison")
        ax[2,0].set_ylabel("Magnitude Diff.")

        fig.tight_layout()
        # fig.text(0.25, 0, s=f"Riemann : {riemann}, Reconst. : {reconst}\nRef Riemann : {ref_riemann}, Ref Reconst. : {ref_reconst}", fontsize=16)
        fig.suptitle(f"t = {round(t_array[ind],4)}\nRiemann : {riemann}, Reconst. : {reconst} --> Ref Riemann : {ref_riemann}, Ref Reconst. : {ref_reconst}", fontsize=16, y=1.0)

        # print(f"Saving to : {self.runname}_{str(ind).zfill(5)}_diff_plot.png")
        plt.savefig(f"{self.path}/{self.runname}/{self.runname}_{riemann}_{reconst}_{str(ind).zfill(5)}_diff_plot.png",dpi=276)
        plt.close()


    def getFinalSlices(self, ind, ds, ref_ds, riemann, reconst, ref_riemann, ref_reconst):

        start = time.time()
        for field in ds.field_list:   

            title = f"{self.runname}_{riemann}_{reconst}_{field[-1]}_{str(ind).zfill(5)}"
            reftitle = f"{self.runname}_{ref_riemann}_{ref_reconst}_{field[-1]}_{str(ind).zfill(5)}"

            p = yt.SlicePlot(ds, "z", fields=field)
            p.annotate_timestamp(time_format="t = {time:.3f}", text_args={'color':'white',
                                'horizontalalignment':'center', 'verticalalignment':'top', 
                                "size":20})

            p.annotate_title(title)
            p.save(name=f"{self.path}/{self.runname}/{title}.png")

            q = yt.SlicePlot(ref_ds, "z", fields=field)
            q.annotate_timestamp(time_format="t = {time:.3f}", text_args={'color':'white',
                                'horizontalalignment':'center', 'verticalalignment':'top', 
                                "size":20})

            q.annotate_title(reftitle)
            q.save(name=f"{self.path}/{self.runname}/{reftitle}.png")

        end = time.time()
        with open(f"{self.path}/{self.runname}/SLICES_FINISHED.txt", "w+") as f:

                f.write(f"Finished with slices in {end - start} sec. \n")
                f.close()

    def checkFinalSlices(self):

        return os.path.isfile(f"{self.path}/{self.savename}/SLICES_FINISHED.txt")

    def getNorms(self, l1_lists, l2_lists, linf_lists, t_lists, riemanns, reconsts, field_list):

        start = time.time()

        cols = plt.cm.cool(np.linspace(0,1,len(l1_lists)))

        for field in range(len(field_list)):

            fig, ax = plt.subplots(3, 1, sharex=True, figsize=(8,10))

            for i in range(len(t_lists)):

                ax[0].plot(t_lists[i], l1_lists[i][field,:], color=cols[i], label=f"{riemanns[i]}_{reconsts[i]}")

                ax[1].plot(t_lists[i], l2_lists[i][field,:], color=cols[i], label=f"{riemanns[i]}_{reconsts[i]}")
            
                ax[2].plot(t_lists[i], linf_lists[i][field,:], color=cols[i], label=f"{riemanns[i]}_{reconsts[i]}")

            ax[0].set_ylabel(r"$L_1$", rotation=0, labelpad = 10)

            ax[1].set_ylabel(r"$L_2$", rotation=0, labelpad = 10)

            ax[2].set_ylabel(r"$L_{\infty}$", rotation=0, labelpad = 10)
            ax[2].set_xlabel(r"Time [code time]")

            for p in range(3):
                ax[p].legend()
                ax[p].grid()

            savetitle = f"{self.path}/{self.savename}/{field_list[field][-1]}_norms_vs_time.png"

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

        self.load_params()

        runname = self.runname
        refname = f"REF_{runname}"

        print(self.runname)

        runs  = glob.glob(f"{self.simpath}/{self.runname}_*")
        refrun = glob.glob(f"{self.simpath}/REF_{self.runname}_*")[0]


        l1_lists   = []
        l2_lists   = []
        linf_lists = []
        t_lists    = []
        riemanns   = []
        recons     = []
        field_list = 0

        for run in runs:
            
            print(f"\n[STARTING] Starting analysis for {run}")

            outputs = sorted(glob.glob(f"{run}/*.phdf"))
            ref_out = sorted(glob.glob(f"{refrun}/*.phdf"))

            riemann = outputs[0].split("/")[-2].split("_")[-2]
            reconst = outputs[0].split("/")[-2].split("_")[-1]

            ref_riemann = ref_out[0].split("/")[-2].split("_")[-2]
            ref_reconst = ref_out[0].split("/")[-2].split("_")[-1]

            temp_val = yt.load(outputs[-1])
            field_list = temp_val.field_list

            num_fields = len(field_list)
            num_touts  = len(outputs)

            l1_array   = np.zeros((num_fields, num_touts), dtype=np.float64)
            l2_array   = np.zeros((num_fields, num_touts), dtype=np.float64)
            linf_array = np.zeros((num_fields, num_touts), dtype=np.float64)
            t_array    = np.zeros(num_touts, dtype=np.float64)

            for ind in tqdm(range(len(outputs)), desc=f"Analyzing {run}"):

                ds = yt.load(outputs[ind])
                refds = yt.load(ref_out[ind])

                # self.getNormVals(ind, ds, refds, l1_array, l2_array, linf_array, t_array)

                self.getDifferencePlots(ind, ds, refds, t_array, riemann, reconst, ref_riemann, ref_reconst)

                # if (ind == 0) or (ind == len(outputs)-1):

                #     self.getFinalSlices(ind, ds, refds, riemann, reconst, ref_riemann, ref_reconst)

            l1_lists.append(l1_array)
            l2_lists.append(l2_array)
            linf_lists.append(linf_array)
            t_lists.append(t_array)
            riemanns.append(riemann)
            recons.append(reconst)

        # self.getNorms(l1_lists, l2_lists, linf_lists, t_lists, riemanns, recons, field_list)

        # # Check if slices complete. If not, run getFinalSlices

        # print("Getting slices")
        # slice_comp = self.checkFinalSlices()
        # if slice_comp == True:
        #     print("Already finished getting slices!")
        #     pass
        # else:
        #     print("RUNNING SLICES")
        #     self.getFinalSlices()
        #     print("Finished getting slices")

        # # Check if norms plot complete. If not, run getNorms

        # print("Getting norms")
        # norm_comp = self.checkNorms()
        # if norm_comp == True:
        #     print("Already finished getting norms!")
        #     pass
        # else:
        #     print("RUNNING NORMS")
        #     self.getNorms()
        #     print("Finished getting norms")

        # Check if final errors compelte. If not, run getFinalErrors



        # Generate HTML file/webpage



        return 0




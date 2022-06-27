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
from matplotlib.colors import SymLogNorm

yt.set_log_level(50)
yt.enable_parallelism()

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

    def getDivB(self, ds):

        # Assumes periodic boundary conditions

        ds.add_gradient_fields(("parthenon", "MagneticField1"))
        ds.add_gradient_fields(("parthenon", "MagneticField2"))
        ds.add_gradient_fields(("parthenon", "MagneticField3"))
        ds.force_periodicity()
        
        ad = ds.all_data()
        dBdx = ad[("parthenon","MagneticField1_gradient_x")].to_ndarray().reshape(ds.domain_dimensions)
        dBdy = ad[("parthenon","MagneticField2_gradient_y")].to_ndarray().reshape(ds.domain_dimensions)
        dBdz = ad[("parthenon","MagneticField3_gradient_z")].to_ndarray().reshape(ds.domain_dimensions)
        
        divB = dBdx + dBdy + dBdz
        
        return np.rot90(divB)

    def extractArrayDivB(self, ds, field):

        ad = ds.all_data()
        arr = ad[field]
        arr_np = arr.to_ndarray().reshape(ds.domain_dimensions)
        
        return arr_np[:,:,0]

    def getcax(self, axis):
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        return cax

    def makeDivBPlot(self, ind, ds, riemann, reconst):

        p_rho = yt.ProjectionPlot(ds, "z", "Density")
        field="Density"
        p_rho.set_cmap(field,"viridis")
        p_rho.set_log(field, log=True)

        p_u   = yt.ProjectionPlot(ds, "z", "Velocity1")
        field="Velocity1"
        p_u.set_cmap(field,"viridis")
        p_u.set_log(field, log=True)

        p_v   = yt.ProjectionPlot(ds, "z", "Velocity2")
        field="Velocity2"
        p_v.set_cmap(field,"viridis")
        p_v.set_log(field, log=True)

        p_Bx  = yt.ProjectionPlot(ds, "z", "MagneticField1")
        field="MagneticField1"
        p_Bx.set_cmap(field,"jet")
        p_Bx.set_log(field, log=True)

        p_By  = yt.ProjectionPlot(ds, "z", "MagneticField2")
        field="MagneticField2"
        p_By.set_cmap(field,"jet")
        p_By.set_log(field, log=True)

        rho  = self.extractArrayDivB(ds, "Density")
        u    = self.extractArrayDivB(ds, "Velocity1")
        v    = self.extractArrayDivB(ds, "Velocity2")
        Bx   = self.extractArrayDivB(ds, "MagneticField1")
        By   = self.extractArrayDivB(ds, "MagneticField2")
        divB = self.getDivB(ds)[:,:,0]
        
        fig, ax = plt.subplots(2, 3, figsize=(16,10))

        im1 = ax[0,0].imshow(rho, cmap="viridis")
        cax1 = self.getcax(ax[0,0])
        plt.colorbar(im1, cax=cax1)
        ax[0,0].set_title(r"$\rho$",fontsize=18, pad=10)
        ax[0,0].axis("off")

        im2 = ax[0,1].imshow(u, cmap="viridis")
        cax2 = self.getcax(ax[0,1])
        plt.colorbar(im2, cax=cax2)
        ax[0,1].set_title(r"$u$",fontsize=18, pad=10)
        ax[0,1].axis("off")

        im3 = ax[0,2].imshow(v, cmap="viridis")
        cax3 = self.getcax(ax[0,2])
        plt.colorbar(im3, cax=cax3)
        ax[0,2].set_title(r"$v$",fontsize=18, pad=10)
        ax[0,2].axis("off")

        im4 = ax[1,0].imshow(Bx, cmap="jet")
        cax4 = self.getcax(ax[1,0])
        plt.colorbar(im4, cax=cax4)
        ax[1,0].set_title(r"$B_x$",fontsize=18, pad=10)
        ax[1,0].axis("off")

        im5 = ax[1,1].imshow(By, cmap="jet", norm=SymLogNorm(linthresh=0.01,vmin=By.min(), vmax=By.max()))
        cax5 = self.getcax(ax[1,1])
        plt.colorbar(im5, cax=cax5)
        ax[1,1].set_title(r"$B_y$",fontsize=18, pad=10)
        ax[1,1].axis("off")

        im6 = ax[1,2].imshow(divB, cmap="jet", norm=SymLogNorm(linthresh=0.1,vmin=divB.min(), vmax=divB.max()))
        cax6 = self.getcax(ax[1,2])
        plt.colorbar(im6, cax=cax6)
        ax[1,2].set_title(r"$\nabla \cdot B$",fontsize=18, pad=10)
        ax[1,2].axis("off")

        fig.suptitle(f"t = {round(float(ds.current_time.value),4)}\nRiemann : {riemann}, Reconst. : {reconst}", fontsize=16, y=1.0)
        
        plt.savefig(f"{self.path}/{self.runname}/{self.runname}_{riemann}_{reconst}_{str(ind).zfill(5)}_divB_plot.png",dpi=276)
        plt.close()

    def getNormVals(self, ind, ds, refds, l1_array, l2_array, linf_array, t_array, normalized=True):

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

            if normalized:

                e_ij = np.abs(np.abs(ref_data[np.where(ref_data != 0.0)]) - np.abs(data[np.where(ref_data != 0.0)])) / np.abs(ref_data[np.where(ref_data != 0.0)])

            else:

                e_ij = np.abs(np.abs(ref_data) - np.abs(data))

            l1 = 1/N * (e_ij.sum())
            l2 = np.sqrt(1/N * ((e_ij**2).sum()))

            try:
                # For initial velocities, all values are zero --> numpy doesn't like this so set linf to 0
                linf = np.amax(e_ij)
            except ValueError:
                linf = 0.0

            l1_array[j, ind] = l1
            l2_array[j, ind] = l2
            linf_array[j, ind] = linf
            t_array[ind] = t

    def getDifferencePlots(self, ind, ds, refds, t_array, riemann, reconst, ref_riemann, ref_reconst):

        ad = ds.all_data()
        ref_ad = refds.all_data()

        num_fields = len(ds.field_list)
        

        MHD  = 0
        TWOD = 1

        for field in ds.field_list:
            # Only MHD runs have Magnetic fields
            if field[-1] == "MagneticPhi":
                MHD = 1

        if ds.domain_dimensions[-1] > 1:
            TWOD = 0

        if TWOD == 1 and MHD == 1:
            num_fields -= 3
        elif TWOD ==1 and MHD == 0:
            num_fields -= 1
        else:
            print("NEED TO IMPLEMENT 3D later")

        if MHD == 0:

            fig, ax = plt.subplots(3, num_fields, figsize=(15,9))

            new_fields = []

            for field in ds.field_list:
                if (field[-1] != "Velocity3"):
                    new_fields.append(field)

            for (j,field) in enumerate(new_fields):
                
                if (TWOD == 1) and (field[-1] != "Velocity3"):
                    arr = ad[field]
                    ref_arr = ref_ad[field]

                    data = arr.to_ndarray().reshape(refds.domain_dimensions)
                    ref_data = ref_arr.to_ndarray().reshape(refds.domain_dimensions)

                    diff = ref_data - data

                    vmax = max(np.amax(data), np.amax(ref_data))
                    vmin = min(np.amin(data), np.amin(ref_data))

                    im1 = ax[0,j].imshow(ref_data[:,:,0], vmax=vmax, vmin=vmin)
                    im2 = ax[1,j].imshow(data[:,:,0], vmax=vmax, vmin=vmin)
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

        else:

            fig, ax = plt.subplots(3, num_fields, figsize=(20,11))

            new_fields = []

            for field in ds.field_list:
                if (field[-1] != "Velocity3") and (field[-1] != "MagneticField3") and (field[-1] != "MagneticPhi"):
                    new_fields.append(field)

            for (j,field) in enumerate(new_fields):
                # print(j, field, num_fields)
                
                if (TWOD == 1) and ((field[-1] != "Velocity3") and (field[-1] != "MagneticField3") and (field[-1] != "MagneticPhi")):
                    arr = ad[field]
                    ref_arr = ref_ad[field]

                    data = arr.to_ndarray().reshape(refds.domain_dimensions)
                    ref_data = ref_arr.to_ndarray().reshape(refds.domain_dimensions)

                    diff = ref_data - data

                    vmax = max(np.amax(data), np.amax(ref_data))
                    vmin = min(np.amin(data), np.amin(ref_data))

                    im1 = ax[0,j].imshow(ref_data[:,:,0], vmax=vmax, vmin=vmin)
                    im2 = ax[1,j].imshow(data[:,:,0], vmax=vmax, vmin=vmin)
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

        wholestart = time.time()

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

            # Make plot with divB?
            divB = True

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

                if divB:
                    self.makeDivBPlot(ind, ds, riemann, reconst)

                self.getNormVals(ind, ds, refds, l1_array, l2_array, linf_array, t_array, normalized=True)

                self.getDifferencePlots(ind, ds, refds, t_array, riemann, reconst, ref_riemann, ref_reconst)

                if (ind == 0) or (ind == len(outputs)-1):

                    self.getFinalSlices(ind, ds, refds, riemann, reconst, ref_riemann, ref_reconst)

            l1_lists.append(l1_array)
            l2_lists.append(l2_array)
            linf_lists.append(linf_array)
            t_lists.append(t_array)
            riemanns.append(riemann)
            recons.append(reconst)

        self.getNorms(l1_lists, l2_lists, linf_lists, t_lists, riemanns, recons, field_list)

        wholeend = time.time()

        print(f"[FINISHED] Finished analysis in {wholeend-wholestart} sec. ({(wholeend-wholestart)/60} min.)")
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




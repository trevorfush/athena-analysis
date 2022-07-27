from mimetypes import init
from weakref import ref
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
import multiprocessing as mp 
import ctypes
import h5py


yt.set_log_level(50)
yt.enable_parallelism()

class Analysis():

    def __init__(self, params, simloc, outfile, saveim):


        # RUNNAME IS THE NAME IN THE CONFIG FILE!!!!
        self.runname   = params["runname"]

        # PROBLEM NAME IS THE NAME OF THE CORRESPONDING INPUT FILE IN ATHENA-PK
        self.inputfile = params["problem_name"]

        # PROBLEM ID IS THE PROBLEM ID FROM THE INPUT FILE FOR ATHENA-PK
        self.problem_id = params["problem_id"]

        self.simloc = simloc

        self.recon     = params["reconstruction"]
        self.riemann   = params["riemann"]

        self.reference = params["reference_sol"]

        self.path      = "results"

        self.kathenasuffix  = "athdf"
        self.athenapksuffix = "phdf"

        self.simpath   = "data"
        self.slice_ind = 0

        self.saveim = saveim

        # Creating field mapping since fields have different names
        self.kathena_field_dict = {
            "Density"        : ('athena_pp','rho'),
            "Velocity1"      : ('athena_pp','vel1'),
            "Velocity2"      : ('athena_pp','vel2'),
            "Velocity3"      : ('athena_pp','vel3'),
            "MagneticField1" : ('athena_pp','Bcc1'),
            "MagneticField2" : ('athena_pp','Bcc2'),
            "MagneticField3" : ('athena_pp','Bcc3'),
            "Pressure"       : ('athena_pp','press'),
            "grad_x"         : ('athena_pp','Bcc1_gradient_x'),
            "grad_y"         : ('athena_pp','Bcc2_gradient_y'),
            "grad_z"         : ('athena_pp','Bcc3_gradient_z')
        }

        self.kathena_field_dict_short = {
            "Density"        : "rho",
            "Velocity1"      : "vel1",
            "Velocity2"      : "vel2",
            "Velocity3"      : "vel3",
            "MagneticField1" : "Bcc1",
            "MagneticField2" : "Bcc2",
            "MagneticField3" : "Bcc3",
            "Pressure"       : "press",
            "grad_x"         : "Bcc1_gradient_x",
            "grad_y"         : "Bcc2_gradient_y",
            "grad_z"         : "Bcc3_gradient_z"
        }

        # Creating field mapping since fields have different names
        self.athenapk_field_dict = {
            "Density"        : ('parthenon','Density'),
            "Velocity1"      : ('parthenon','Velocity1'),
            "Velocity2"      : ('parthenon','Velocity2'),
            "Velocity3"      : ('parthenon','Velocity3'),
            "MagneticField1" : ('parthenon','MagneticField1'),
            "MagneticField2" : ('parthenon','MagneticField2'),
            "MagneticField3" : ('parthenon','MagneticField3'),
            "Pressure"       : ('parthenon','Pressure'),
            "grad_x"         : ('parthenon','MagneticField1_gradient_x'),
            "grad_y"         : ('parthenon','MagneticField2_gradient_y'),
            "grad_z"         : ('parthenon','MagneticField3_gradient_z')
        }

        self.athenapk_field_dict_short = {
            "Density"        : "Density",
            "Velocity1"      : "Velocity1",
            "Velocity2"      : "Velocity2",
            "Velocity3"      : "Velocity3",
            "MagneticField1" : "MagneticField1",
            "MagneticField2" : "MagneticField2",
            "MagneticField3" : "MagneticField3",
            "Pressure"       : "Pressure",
            "grad_x"         : "MagneticField1_gradient_x",
            "grad_y"         : "MagneticField2_gradient_y",
            "grad_z"         : "MagneticField3_gradient_z"
        }

        self.athenapk  = "../athenapk"

        if self.reference == True:
            self.savename = f"REF_{self.runname}"
        else:
            self.savename = self.runname

        
        self.plotsavepreamble = f"{self.runname}_{self.riemann}_{self.recon}"
        
        self.simparams = 0

        self.outfile = outfile

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
            current_progress[2] = max(glob.glob(f"{self.simpath}/{self.savename}/*.?????.{self.filesuffix}"))
        if (fileexist == False and direxist == False):
            current_progress[1] = False
            current_progress[2] = None

        return current_progress

    def getDivB(self, ds, field_dict):

        # Assumes periodic boundary conditions

        ds.add_gradient_fields(field_dict["MagneticField1"])
        ds.add_gradient_fields(field_dict["MagneticField2"])
        ds.add_gradient_fields(field_dict["MagneticField3"])
        ds.force_periodicity()
        
        ad = ds.all_data()
        dBdx = ad[field_dict["grad_x"]].to_ndarray().reshape(ds.domain_dimensions)
        dBdy = ad[field_dict["grad_y"]].to_ndarray().reshape(ds.domain_dimensions)
        dBdz = ad[field_dict["grad_z"]].to_ndarray().reshape(ds.domain_dimensions)
        
        divB = dBdx + dBdy + dBdz
        
        return np.rot90(divB)

    def extractArrayDivB(self, ds, field):

        ad = ds.all_data()
        arr = ad[field]
        arr_np = arr.to_ndarray().reshape(ds.domain_dimensions)
        
        return np.rot90(arr_np[:,:,self.slice_ind])

    def getcax(self, axis):
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        return cax

    def makeDivBPlot(self, id, code, ind, ds, riemann, reconst, field_dict_list, field_dict_short_list):

        

        field_dict       = field_dict_list[id]
        field_dict_short = field_dict_short_list[id]

        if self.runname == "turbulence":
            c = [0.5, 0.5, 1/3]
        else:
            c = "c"

        # print(f"CODE for ind {ind} = {code}")
        # print(f"FIELD DICT short : {field_dict}")

        if self.saveim == True:
            p_rho = yt.SlicePlot(ds, "z", field_dict_short["Density"], center=c)
            field=field_dict_short["Density"]
            p_rho.set_cmap(field,"viridis")
            p_rho.set_log(field, log=True)

            p_u   = yt.SlicePlot(ds, "z", field_dict_short["Velocity1"], center=c)
            field=field_dict_short["Velocity1"]
            p_u.set_cmap(field,"viridis")
            p_u.set_log(field, log=True)

            p_v   = yt.SlicePlot(ds, "z", field_dict_short["Velocity2"], center=c)
            field=field_dict_short["Velocity2"]
            p_v.set_cmap(field,"viridis")
            p_v.set_log(field, log=True)

            p_Bx  = yt.SlicePlot(ds, "z", field_dict_short["MagneticField1"], center=c)
            field=field_dict_short["MagneticField1"]
            p_Bx.set_cmap(field,"jet")
            p_Bx.set_log(field, log=True)

            p_By  = yt.SlicePlot(ds, "z", field_dict_short["MagneticField2"], center=c)
            field=field_dict_short["MagneticField2"]
            p_By.set_cmap(field,"jet")
            p_By.set_log(field, log=True)

            # rho  = self.extractArrayDivB(ds, self.field_dict_short["Density"])
            # u    = self.extractArrayDivB(ds, self.field_dict_short["Velocity1"])
            # v    = self.extractArrayDivB(ds, self.field_dict_short["Velocity2"])
            # Bx   = self.extractArrayDivB(ds, self.field_dict_short["MagneticField1"])
            # By   = self.extractArrayDivB(ds, self.field_dict_short["MagneticField2"])
            # divB = self.getDivB(ds)[:,:,self.slice_ind]

            rho  = np.array(p_rho.frb[field_dict["Density"]])
            u    = np.array(p_rho.frb[field_dict["Velocity1"]])
            v    = np.array(p_rho.frb[field_dict["Velocity2"]])
            Bx   = np.array(p_rho.frb[field_dict["MagneticField1"]])
            By   = np.array(p_rho.frb[field_dict["MagneticField2"]])

        temp_divB = self.getDivB(ds, field_dict)
        divB = temp_divB[:,:,self.slice_ind]

        mean_divB = temp_divB.mean()
        max_divB = np.abs(temp_divB).max()
        
        if self.saveim == True:
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

            fig.suptitle(f"t = {round(float(ds.current_time.value),4)}\nRiemann : {riemann}, Reconst. : {reconst}, Code : {code}", fontsize=16, y=1.0)
            
            plt.savefig(f"{self.path}/{self.runname}/{self.runname}_{code}_{riemann}_{reconst}_{str(ind).zfill(5)}_divB_plot.png",dpi=276)
            plt.close()

        return max_divB, mean_divB

    def getNormVals(self, id, code, ind, ds, refds, l1_array, l2_array, linf_array, t_array, field_dict_list, field_dict_short_list, normalized=True):

        # Getting reference riemann and reconstruction
        # reference_strs = glob.glob(f"{self.simpath}/REF_{self.runname}_*")
        # reference_str  = reference_strs[0]

        # ref_riemann    = reference_str.split("_")[-2]
        # ref_recon      = reference_str.split("_")[-1]

        # ds     = yt.load(output)
        # ref_ds = yt.load(f"{self.simpath}/REF_{self.runname}_{ref_riemann}_{ref_recon}/{output.split('/')[-1]}")
        field_dict       = field_dict_list[id]
        field_dict_short = field_dict_short_list[id]

        t = np.float64(ds.current_time.value)

        cg  = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
        refcg  = refds.covering_grid(level=0, left_edge=refds.domain_left_edge, dims=refds.domain_dimensions)

        new_fields = []

        for field in ds.field_list:
            if (field[-1] != "MagneticPhi"):
                new_fields.append(field)

        for (j, field) in enumerate(new_fields):

            data = cg[field].d

            if (code == "athenapk" and self.refcode == "athenapk") or (code == "kathena" and self.refcode == "kathena"):
                refarr = refcg[field]
            if code == "athenapk" and self.refcode == "kathena":
                pass
            if code == "kathena" and self.refcode == "athenapk":
                if field == self.kathena_field_dict["Density"]:
                    refarr = refcg[self.athenapk_field_dict["Density"]]
                elif field == self.kathena_field_dict["Velocity1"]:
                    refarr = refcg[self.athenapk_field_dict["Velocity1"]]
                elif field == self.kathena_field_dict["Velocity2"]:
                    refarr = refcg[self.athenapk_field_dict["Velocity2"]]
                elif field == self.kathena_field_dict["Velocity3"]:
                    refarr = refcg[self.athenapk_field_dict["Velocity3"]]
                elif field == self.kathena_field_dict["MagneticField1"]:
                    refarr = refcg[self.athenapk_field_dict["MagneticField1"]]
                elif field == self.kathena_field_dict["MagneticField2"]:
                    refarr = refcg[self.athenapk_field_dict["MagneticField2"]]
                elif field == self.kathena_field_dict["MagneticField3"]:
                    refarr = refcg[self.athenapk_field_dict["MagneticField3"]]
                elif field == self.kathena_field_dict["Pressure"]:
                    refarr = refcg[self.athenapk_field_dict["Pressure"]]
                else:
                    print("Whoops, check here for issue")

            ref_data = refarr.d

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
                # For some initial fields, all values are zero --> numpy doesn't like this so set linf to 0
                linf = np.amax(e_ij)
            except ValueError:
                linf = 0.

            l1_array[j, ind] = l1
            l2_array[j, ind] = l2
            linf_array[j, ind] = linf
            t_array[ind] = t

    def getDifferencePlots(self, id, code, ind, ds, refds, t_array, riemann, reconst, ref_riemann, ref_reconst, field_dict_list, field_dict_short_list):

        field_dict       = field_dict_list[id]
        field_dict_short = field_dict_short_list[id]

        # ad = ds.all_data()
        # ref_ad = refds.all_data()

        num_fields = len(ds.field_list)
        

        MHD  = 0
        TWOD = 1

        for field in ds.field_list:
            # Only MHD runs have Magnetic fields
            if field[-1] == field_dict_short["MagneticField1"] or self.runname == "turbulence":
                MHD = 1

        if ds.domain_dimensions[-1] > 1 or self.runname == "turbulence":
            TWOD = 0

        if TWOD == 1 and MHD == 1:
            num_fields -= 3
        elif TWOD ==1 and MHD == 0:
            num_fields -= 1
        else:
            print("Make sure you're running turbulence")

        if MHD == 0 and TWOD == 1:

            # fig, ax = plt.subplots(3, num_fields, figsize=(15,9))

            new_fields = []

            # for field in ds.field_list:
            #     if (field[-1] != field_dict_short["Velocity3"]):
            #         new_fields.append(field)

            # num_fields == len(new_fields)

            # for (j,field) in enumerate(new_fields):
                
            #     if (TWOD == 1) and (field[-1] != field_dict_short["Velocity3"]):
            #         arr = ad[field]

            #         if (code == "athenapk" and self.refcode == "athenapk") or (code == "kathena" and self.refcode == "kathena"):
            #             refarr = ref_ad[field]
            #         if code == "athenapk" and self.refcode == "kathena":
            #             pass
            #         if code == "kathena" and self.refcode == "athenapk":
            #             if field == self.kathena_field_dict["Density"]:
            #                 refarr = ref_ad[self.athenapk_field_dict["Density"]]
            #             elif field == self.kathena_field_dict["Velocity1"]:
            #                 refarr = ref_ad[self.athenapk_field_dict["Velocity1"]]
            #             elif field == self.kathena_field_dict["Velocity2"]:
            #                 refarr = ref_ad[self.athenapk_field_dict["Velocity2"]]
            #             elif field == self.kathena_field_dict["Velocity3"]:
            #                 refarr = ref_ad[self.athenapk_field_dict["Velocity3"]]
            #             elif field == self.kathena_field_dict["MagneticField1"]:
            #                 refarr = ref_ad[self.athenapk_field_dict["MagneticField1"]]
            #             elif field == self.kathena_field_dict["MagneticField2"]:
            #                 refarr = ref_ad[self.athenapk_field_dict["MagneticField2"]]
            #             elif field == self.kathena_field_dict["MagneticField3"]:
            #                 refarr = ref_ad[self.athenapk_field_dict["MagneticField3"]]
            #             elif field == self.kathena_field_dict["Pressure"]:
            #                 refarr = ref_ad[self.athenapk_field_dict["Pressure"]]
            #             else:
            #                 print("Whoops, check here for issue")

            #         # ref_arr = ref_ad[field]
            #         if code == "athenapk":
            #             data = arr.to_ndarray().reshape(refds.domain_dimensions)
            #         elif code == "kathena":
            #             data_p = yt.SlicePlot(ds, "z", field, center=c)
            #             data = np.array(data_p.frb[field])

            #         ref_data = refarr.to_ndarray().reshape(refds.domain_dimensions)

            #         diff = ref_data - data

            #         vmax = max(np.amax(data), np.amax(ref_data))
            #         vmin = min(np.amin(data), np.amin(ref_data))

            #         im1 = ax[0,j].imshow(ref_data[:,:,0], vmax=vmax, vmin=vmin)
            #         im2 = ax[1,j].imshow(data[:,:,0], vmax=vmax, vmin=vmin)
            #         im3 = ax[2,j].imshow(diff[:,:,0], cmap="jet")

            #         divider1 = make_axes_locatable(ax[0,j])
            #         cax1 = divider1.append_axes("right",size="5%",pad=0.05)
            #         cbar1 = plt.colorbar(im1, cax=cax1)

            #         divider2 = make_axes_locatable(ax[1,j])
            #         cax2 = divider2.append_axes("right",size="5%",pad=0.05)
            #         cbar2 = plt.colorbar(im2, cax=cax2)

            #         divider3 = make_axes_locatable(ax[2,j])
            #         cax3 = divider3.append_axes("right",size="5%",pad=0.05)
            #         cbar3 = plt.colorbar(im3, cax=cax3)

            #         ax[0,j].set_title(field[-1])

        elif MHD == 1 and TWOD == 1:

            if self.runname == "turbulence":
                c = [0.5, 0.5, 1/3]
            else:
                c = "c"

            new_fields = []

            for field in ds.field_list:
                if (field[-1] != field_dict_short["Velocity3"]) and (field[-1] != field_dict_short["MagneticField3"]) and (field[-1] != "MagneticPhi"):
                    new_fields.append(field)

            num_fields = len(new_fields)

            fig, ax = plt.subplots(3, num_fields, figsize=(20,11))
            # if id == 5:
            #     print(new_fields)
            for (j,field) in enumerate(new_fields):
                # if id == 5:
                #     print(j, field, num_fields)
            
                # arr = ad[field]

                if (code == "athenapk" and self.refcode == "athenapk") or (code == "kathena" and self.refcode == "kathena"):
                    reffield = field
                if code == "athenapk" and self.refcode == "kathena":
                    pass
                if code == "kathena" and self.refcode == "athenapk":
                    if field == self.kathena_field_dict["Density"]:
                        reffield = self.athenapk_field_dict["Density"]
                    elif field == self.kathena_field_dict["Velocity1"]:
                        reffield = self.athenapk_field_dict["Velocity1"]
                    elif field == self.kathena_field_dict["Velocity2"]:
                        reffield = self.athenapk_field_dict["Velocity2"]
                    elif field == self.kathena_field_dict["Velocity3"]:
                        reffield = self.athenapk_field_dict["Velocity3"]
                    elif field == self.kathena_field_dict["MagneticField1"]:
                        reffield = self.athenapk_field_dict["MagneticField1"]
                    elif field == self.kathena_field_dict["MagneticField2"]:
                        reffield = self.athenapk_field_dict["MagneticField2"]
                    elif field == self.kathena_field_dict["MagneticField3"]:
                        reffield = self.athenapk_field_dict["MagneticField3"]
                    elif field == self.kathena_field_dict["Pressure"]:
                        reffield = self.athenapk_field_dict["Pressure"]
                    else:
                        print("Whoops, check here for issue")

                # ref_arr = ref_ad[field]
                data_p = yt.SlicePlot(ds, "z", field, center=c)
                data = np.array(data_p.frb[field])

                refdata_p = yt.SlicePlot(refds, "z", reffield, center=c)
                ref_data = np.array(refdata_p.frb[reffield])

               #  data = arr.to_ndarray().reshape(refds.domain_dimensions)
               # ref_data = refarr.to_ndarray().reshape(refds.domain_dimensions)

                diff = ref_data - data

                vmax = max(np.amax(data), np.amax(ref_data))
                vmin = min(np.amin(data), np.amin(ref_data))

                im1 = ax[0,j].imshow(ref_data[:,:], vmax=vmax, vmin=vmin)
                im2 = ax[1,j].imshow(data[:,:], vmax=vmax, vmin=vmin)
                im3 = ax[2,j].imshow(diff[:,:], cmap="jet")

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

        # This should apply to only the turbulence runs
        elif MHD == 1 and TWOD == 0:

            if self.runname == "turbulence":
                c = [0.5, 0.5, 1/3]
            else:
                c = "c"

            new_fields = []

            for field in ds.field_list:
                new_fields.append(field)

            num_fields = len(new_fields)

            fig, ax = plt.subplots(3, num_fields, figsize=(20,11))

            for (j,field) in enumerate(new_fields):
                
                # arr = ad[field]
                # ref_arr = ref_ad[field]

                # data = arr.to_ndarray().reshape(refds.domain_dimensions)
                # ref_data = ref_arr.to_ndarray().reshape(refds.domain_dimensions)

                # diff = ref_data - data

                if (code == "athenapk" and self.refcode == "athenapk") or (code == "kathena" and self.refcode == "kathena"):
                    reffield = field
                if code == "athenapk" and self.refcode == "kathena":
                    pass
                if code == "kathena" and self.refcode == "athenapk":
                    if field == self.kathena_field_dict["Density"]:
                        reffield = self.athenapk_field_dict["Density"]
                    elif field == self.kathena_field_dict["Velocity1"]:
                        reffield = self.athenapk_field_dict["Velocity1"]
                    elif field == self.kathena_field_dict["Velocity2"]:
                        reffield = self.athenapk_field_dict["Velocity2"]
                    elif field == self.kathena_field_dict["Velocity3"]:
                        reffield = self.athenapk_field_dict["Velocity3"]
                    elif field == self.kathena_field_dict["MagneticField1"]:
                        reffield = self.athenapk_field_dict["MagneticField1"]
                    elif field == self.kathena_field_dict["MagneticField2"]:
                        reffield = self.athenapk_field_dict["MagneticField2"]
                    elif field == self.kathena_field_dict["MagneticField3"]:
                        reffield = self.athenapk_field_dict["MagneticField3"]
                    elif field == self.kathena_field_dict["Pressure"]:
                        reffield = self.athenapk_field_dict["Pressure"]
                    else:
                        print("Whoops, check here for issue")

                data_p = yt.SlicePlot(ds, "z", field, center=c)
                ref_data_p = yt.SlicePlot(refds, "z", reffield, center=c)

                data = np.array(data_p.frb[field])
                ref_data = np.array(ref_data_p.frb[reffield])

                diff = data - ref_data

                vmax = max(np.amax(data), np.amax(ref_data))
                vmin = min(np.amin(data), np.amin(ref_data))

                im1 = ax[0,j].imshow(ref_data, vmax=vmax, vmin=vmin)
                im2 = ax[1,j].imshow(data, vmax=vmax, vmin=vmin)
                im3 = ax[2,j].imshow(diff, cmap="jet")

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
        fig.suptitle(f"t = {round(t_array[ind],4)}\nRiemann : {riemann}, Reconst. : {reconst} Code : {code} --> Ref Riemann : {ref_riemann}, Ref Reconst. : {ref_reconst}, Ref code : {self.refcode}", fontsize=16, y=1.0)

        # print(f"Saving to : {self.runname}_{str(ind).zfill(5)}_diff_plot.png")
        plt.savefig(f"{self.path}/{self.runname}/{self.runname}_{code}_{riemann}_{reconst}_{str(ind).zfill(5)}_diff_plot.png",dpi=276)
        plt.close()


    def getFinalSlices(self, id, code, ind, ds, ref_ds, riemann, reconst, ref_riemann, ref_reconst, field_dict_list, field_dict_short_list):


        field_dict       = field_dict_list[id]
        field_dict_short = field_dict_short_list[id]

        if self.runname == "turbulence":
            c = [0.5, 0.5, 1/3]
        else:
            c = "c"

        start = time.time()
        for field in ds.field_list:  

            if (code == "athenapk" and self.refcode == "athenapk") or (code == "kathena" and self.refcode == "kathena"):
                reffield = field
            if code == "athenapk" and self.refcode == "kathena":
                pass
            if code == "kathena" and self.refcode == "athenapk":
                if field == self.kathena_field_dict["Density"]:
                    reffield = self.athenapk_field_dict["Density"]
                elif field == self.kathena_field_dict["Velocity1"]:
                    reffield = self.athenapk_field_dict["Velocity1"]
                elif field == self.kathena_field_dict["Velocity2"]:
                    reffield = self.athenapk_field_dict["Velocity2"]
                elif field == self.kathena_field_dict["Velocity3"]:
                    reffield = self.athenapk_field_dict["Velocity3"]
                elif field == self.kathena_field_dict["MagneticField1"]:
                    reffield = self.athenapk_field_dict["MagneticField1"]
                elif field == self.kathena_field_dict["MagneticField2"]:
                    reffield = self.athenapk_field_dict["MagneticField2"]
                elif field == self.kathena_field_dict["MagneticField3"]:
                    reffield = self.athenapk_field_dict["MagneticField3"]
                elif field == self.kathena_field_dict["Pressure"]:
                    reffield = self.athenapk_field_dict["Pressure"]
                else:
                        print("Whoops, check here for issue")

            title = f"{self.runname}_{code}_{riemann}_{reconst}_{field[-1]}_{str(ind).zfill(5)}"
            reftitle = f"{self.runname}_{self.refcode}_{ref_riemann}_{ref_reconst}_{field[-1]}_{str(ind).zfill(5)}"

            p = yt.SlicePlot(ds, "z", fields=field, center=c)
            p.annotate_timestamp(time_format="t = {time:.3f}", text_args={'color':'white',
                                'horizontalalignment':'center', 'verticalalignment':'top', 
                                "size":20})

            p.annotate_title(title)
            p.save(name=f"{self.path}/{self.runname}/{title}.png")

            q = yt.SlicePlot(ref_ds, "z", fields=reffield, center=c)
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

    def getNorms(self, l1_lists, l2_lists, linf_lists, t_lists, riemanns, reconsts, codes, field_list):

        

        start = time.time()
        
        cols = {
            "hlle_dc"    : "black",
            "hlle_plm"   : "lightcoral",
            "hlle_ppm"   : "firebrick",
            "hlle_wenoz" : "sienna",
            "hlld_dc"    : "forestgreen",
            "hlld_plm"   : "purple",
            "hlld_ppm"   : "dodgerblue",
            "hlld_wenoz" : "deeppink"
        }

        linestyles = {
            "athenapk" : "solid",
            "kathena"  : "dotted"
        }

        new_fields = []

        for field in field_list:
            if (field[-1] != "MagneticPhi"):
                new_fields.append(field)

        for field in range(len(new_fields)):

            fig, ax = plt.subplots(3, 1, sharex=True, figsize=(8,10))

            for i in range(len(t_lists)):

                ax[0].plot(t_lists[i][:], l1_lists[i][field,:], 
                           color = cols[f"{riemanns[i]}_{reconsts[i]}"], 
                           label=f"{codes[i]}_{riemanns[i]}_{reconsts[i]}",
                           linestyle=linestyles[f"{codes[i]}"])

                ax[1].plot(t_lists[i][:], l2_lists[i][field,:], 
                           color = cols[f"{riemanns[i]}_{reconsts[i]}"], 
                           label=f"{codes[i]}_{riemanns[i]}_{reconsts[i]}",
                           linestyle=linestyles[f"{codes[i]}"])
            
                ax[2].plot(t_lists[i][:], linf_lists[i][field,:], 
                           color = cols[f"{riemanns[i]}_{reconsts[i]}"], 
                           label=f"{codes[i]}_{riemanns[i]}_{reconsts[i]}",
                           linestyle=linestyles[f"{codes[i]}"])

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

    def getDivBPlot(self, max_divBs, mean_divBs, ts, recons, riemanns, codes):

        cols = {
            "hlle_dc"    : "black",
            "hlle_plm"   : "lightcoral",
            "hlle_ppm"   : "firebrick",
            "hlle_wenoz" : "sienna",
            "hlld_dc"    : "forestgreen",
            "hlld_plm"   : "purple",
            "hlld_ppm"   : "dodgerblue",
            "hlld_wenoz" : "deeppink"
        }

        linestyles = {
            "athenapk" : "solid",
            "kathena"  : "dotted"
        }

        fig, ax = plt.subplots(2,1, sharex=True, figsize=(8,10))

        for i in range(len(ts)):

            ax[0].plot(ts[i], max_divBs[i], color = cols[f"{riemanns[i]}_{recons[i]}"],
                       label=f"{codes[i]}_{riemanns[i]}_{recons[i]}", linestyle=linestyles[f"{codes[i]}"], linewidth=2)

            ax[1].plot(ts[i], mean_divBs[i], color = cols[f"{riemanns[i]}_{recons[i]}"],
                       label=f"{codes[i]}_{riemanns[i]}_{recons[i]}", linestyle=linestyles[f"{codes[i]}"], linewidth=2)

        ax[0].set_ylabel(r"max($\nabla \cdot B$)")
        ax[1].set_ylabel(r"$\langle \nabla \cdot B \rangle$")
        ax[1].set_xlabel("Time [code time]")

        for p in range(2):
            ax[p].legend()
            ax[p].grid()

        savetitle = f"{self.path}/{self.savename}/divB_vs_time.png"

        plt.savefig(savetitle)

        plt.close()


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

    def writeDataToOutfile(self, runs, l1_arr, l2_arr, linf_arr, max_divB_arr, mean_divB_arr, t_arr):
        """
        Write norm and divB data to outfile (HDF5 file)
        -----------------------------------------------------
        Data file structure:

        File :

            -> Groups represent each of the different runs, 
               noted by the simulation directory name :

               -> Each group has multiple datasets :
                    - L1 norms vs time
                    - L2 norms vs time
                    - Linf norms vs time
                    - Max |divB| vs time
                    - Mean divB vs time
                    - Array of time values
        """

        h5file = h5py.File(f"{self.path}/{self.runname}/{self.outfile}", "w")

        # Make groups for each run
        for i, run in enumerate(runs):
            
            h5group = h5file.create_group(run)

            # L1 norms
            h5dset_l1  = h5group.create_dataset(name="l1norms", shape=l1_arr[i].shape, dtype=np.float64)
            h5dset_l1[:,:] = l1_arr[i][:,:]

            # L2 norms
            h5dset_l2  = h5group.create_dataset(name="l2norms", shape=l2_arr[i].shape, dtype=np.float64)
            h5dset_l2[:,:] = l2_arr[i][:,:]

            # Linf norms
            h5dset_linf  = h5group.create_dataset(name="linfnorms", shape=linf_arr[i].shape, dtype=np.float64)
            h5dset_linf[:,:] = linf_arr[i][:,:]

            # Max div B
            h5dset_maxdivB = h5group.create_dataset(name="maxdivB", shape=max_divB_arr[i].shape, dtype=np.float64)
            h5dset_maxdivB[:] = max_divB_arr[i][:]

            # Mean div B
            h5dset_meandivB = h5group.create_dataset(name="meandivB", shape=mean_divB_arr[i].shape, dtype=np.float64)
            h5dset_meandivB[:] = mean_divB_arr[i][:]

            # t arrays
            h5dset_ts = h5group.create_dataset(name="time", shape=t_arr[i].shape, dtype=np.float64)
            h5dset_ts[:] = t_arr[i][:]

        h5file.close()

        return 0
        

    def runThread(self, id, codes, divB, run, refrun, l1_lists, l2_lists, linf_lists, max_divB_arr, mean_divB_arr, t_lists, field_dict_list, field_dict_short_list):
        """
        Runs analysis pipeline on a given thread for a given combination of numerical methods

        The norm lists are passed through so that values of norms over time can be kept track of
        across processes. 
        """

        code = codes[id]
        print(code, run)

        if code == "athenapk":
            field_dict       = self.athenapk_field_dict
            field_dict_short = self.athenapk_field_dict_short
        else:
            field_dict       = self.kathena_field_dict
            field_dict_short = self.kathena_field_dict_short
            
        # print(id, code)

        # Gather all outputs for the run and the reference solution
        if code == "athenapk":
            if self.runname == "turbulence":
                outputs = sorted(glob.glob(f"{run}/Turb.prim.?????.phdf"))
                if self.refcode == "athenapk":
                    ref_out = sorted(glob.glob(f"{refrun}/Turb.prim.?????.phdf"))
                else:
                    ref_out = sorted(glob.glob(f"{refrun}/Turb.prim.?????.{self.kathenasuffix}"))
                print ("YUP ITS ATHENAPK")
            else:
                outputs = sorted(glob.glob(f"{run}/*.phdf"))
                if self.refcode == "athenapk":
                    ref_out = sorted(glob.glob(f"{refrun}/*.phdf"))
                else:
                    ref_out = sorted(glob.glob(f"{refrun}/*.{self.kathenasuffix}"))
                print ("YUP ITS ATHENAPK")
        else:
            if self.runname == "turbulence":
                outputs = sorted(glob.glob(f"{run}/Turb.prim.?????.{self.kathenasuffix}"))
                if self.refcode == "kathena":
                    ref_out = sorted(glob.glob(f"{refrun}/Turb.prim.?????.{self.kathenasuffix}"))
                else:
                    ref_out = sorted(glob.glob(f"{refrun}/Turb.prim.?????.phdf"))
            else:
                outputs = sorted(glob.glob(f"{run}/*.{self.kathenasuffix}"))
                if self.refcode == "kathena":
                    ref_out = sorted(glob.glob(f"{refrun}/*.{self.kathenasuffix}"))
                else:
                    ref_out = sorted(glob.glob(f"{refrun}/*.phdf"))


        # Get riemann and reconstruction methods for labelling the plots
        
        # print("HERE", run, outputs)
        # print("HERE", run, sorted(glob.glob(f"{run}/*.phdf"))[0])
        riemann = outputs[0].split("/")[-2].split("_")[-2]
        reconst = outputs[0].split("/")[-2].split("_")[-1]
        # code    = outputs[0].split("/")[-2].split("_")[-3]

        ref_riemann = ref_out[0].split("/")[-2].split("_")[-2]
        ref_reconst = ref_out[0].split("/")[-2].split("_")[-1]
        ref_code    = ref_out[0].split("/")[-2].split("_")[-3]

        # Create temporary empty arrays (full of zeros) to store norm info
        # These values will later be added to the true norm array shared 
        # across processors
        temp_val = yt.load(outputs[-1])
        field_list = temp_val.field_list

        num_fields = len(field_list)

        if code == "athenapk":
            num_fields -= 1

        num_touts  = len(outputs)

        print(f"INFO --> code : {code}, num_fields : {num_fields}, num_touts : {num_touts}")

        l1_array   = np.zeros((num_fields, num_touts), dtype=np.float64)
        l2_array   = np.zeros((num_fields, num_touts), dtype=np.float64)
        linf_array = np.zeros((num_fields, num_touts), dtype=np.float64)
        max_divB_array = np.zeros(num_touts, dtype=np.float64)
        mean_divB_array = np.zeros(num_touts, dtype=np.float64)
        t_array    = np.zeros(num_touts, dtype=np.float64)
        
        # Loop through outputs of simulation and conduct analysis
        for ind in tqdm(range(len(outputs)), desc=f"Analyzing {run}"):
            
            # Load comparison and reference dataset
            ds = yt.load(outputs[ind])
            refds = yt.load(ref_out[ind])

            # if code == "kathena":
            #     print(ds.field_list)

            # If the simulation is MHD, make divB plots
            if divB:
                max_divB, mean_divB = self.makeDivBPlot(id, code, ind, ds, riemann, reconst, field_dict_list, field_dict_short_list)
                max_divB_array[ind] = max_divB
                mean_divB_array[ind] = mean_divB

            # Get the norm values for this output and append to the temporary arrays
            self.getNormVals(id, code, ind, ds, refds, l1_array, l2_array, linf_array, t_array, field_dict_list, field_dict_short_list, normalized=True)

            # Make the difference plot for this output
            if self.saveim == True:
                self.getDifferencePlots(id, code, ind, ds, refds, t_array, riemann, reconst, ref_riemann, ref_reconst, field_dict_list, field_dict_short_list)

            # Make slices of final or initial output 
            if self.saveim == True:
                if (ind == 0) or (ind == len(outputs)-1):
                    self.getFinalSlices(id, code, ind, ds, refds, riemann, reconst, ref_riemann, ref_reconst, field_dict_list, field_dict_short_list)

        # Insert the calculated norm values into the correct location in the shared norm lists
        # print(f"L1 array : {l1_array}\nL1 array type : {type(l1_array)}")
        # print(f"L1 list array : {l1_lists[id]}\nL1 array type : {type(l1_lists[id])}")

        # print(f"L1 lists {id} INSIDE: {l1_lists}")

        l1_lists[id][:,:] = l1_array[:,:]
        l2_lists[id][:,:] = l2_array[:,:]
        linf_lists[id][:,:] = linf_array[:,:]
        t_lists[id][:] = t_array[:]

        max_divB_arr[id][:] = max_divB_array[:]
        mean_divB_arr[id][:] = mean_divB_array[:]


    def getEmptyNormArrays(self, runs):
        """
        Generates empty (or zero filled) data structures to store norm and sim
        info across analysis threads. The shapes are:

        l1, l2, linf     --> Number of sims X Number of sim fields X Number of time outputs
        t                --> Number of sims X Number of time outputs
        riemanns, recons --> Number of sims

        Number of sims denotes number of unique combinations of num. methods, 
        number of sim fields is the number of fields in the datasets, and 
        number of time outputs are the number of outputs for each simulation
        """

        l1 = []
        l2 = []
        linf = []
        max_divB = []
        mean_divB = []
        t = []

        for run in runs:
            tempoutput = sorted(glob.glob(f"{run}/*.{self.kathenasuffix}")) + sorted(glob.glob(f"{run}/*.{self.athenapksuffix}"))
            temp_val = yt.load(tempoutput[0])
            field_list = temp_val.field_list

            code    = tempoutput[0].split("/")[-2].split("_")[-3]

            num_fields = len(field_list)
            num_touts  = len(tempoutput)

            if code == "athenapk":
                num_fields -= 1

            l1_lists   = np.zeros((num_fields, num_touts))
            l2_lists   = np.zeros((num_fields, num_touts))
            linf_lists = np.zeros((num_fields, num_touts))
            max_divB_lists = np.zeros((num_touts))
            mean_divB_lists = np.zeros((num_touts))
            t_lists    = np.zeros((num_touts))

            l1.append(l1_lists)
            l2.append(l2_lists)
            linf.append(linf_lists)
            max_divB.append(max_divB_lists)
            mean_divB.append(mean_divB_lists)
            t.append(t_lists)

        riemanns   = np.zeros(len(runs), dtype=np.dtype("U100"))
        recons     = np.zeros(len(runs), dtype=np.dtype("U100"))
        codes      = np.zeros(len(runs), dtype=np.dtype("U100"))
        integs     = np.zeros(len(runs), dtype=np.dtype("U100"))

        return np.array(l1), np.array(l2), np.array(linf), np.array(max_divB), np.array(mean_divB), np.array(t), np.array(riemanns), np.array(recons), np.array(codes), np.array(integs)


    def Analyze(self):
        """
        Runs analysis on given simulation and tracks progress.
        """

        if os.path.isdir(f"{self.path}/{self.savename}") != True:
            os.mkdir(f"{self.path}/{self.savename}")

        self.load_params()

        runname = self.runname
        refname = f"REF_{runname}"

        print(self.runname)

        pre_runs  = glob.glob(f"{self.simloc}/{self.runname}_*")
        pre_refrun = glob.glob(f"{self.simloc}/REF_{self.runname}_*")

        runs = []

        for run in pre_runs:
            if run.split(".")[-1] != "gz":
                runs.append(run)
        
        for run in pre_refrun:
            if run.split(".")[-1] != "gz":
                refrun = run

        self.refcode = refrun.split("/")[-1].split("_")[-3]

        print(f"Runs = {runs}")

        # Generating empty arrays for norms
        ###############################################################
        if self.runname == "turbulence":
            tempoutput = sorted(glob.glob(f"{runs[0]}/Turb.prim.?????.{self.athenapksuffix}")) + sorted(glob.glob(f"{runs[0]}/Turb.prim.?????.{self.kathenasuffix}"))
        else:
            tempoutput = sorted(glob.glob(f"{runs[0]}/*.{self.athenapksuffix}")) + sorted(glob.glob(f"{runs[0]}/*.{self.kathenasuffix}"))

        temp_val = yt.load(tempoutput[-1])
        field_list = temp_val.field_list

        num_fields = len(field_list)
        num_touts  = len(tempoutput)

        l1_lists, l2_lists, linf_lists, max_divB_lists, mean_divB_lists, t_lists, riemanns, recons, codes, integs = self.getEmptyNormArrays(runs)

        ###############################################################

        # print(f"L1 list size : {l1_lists.shape}")

        norm_array_sizes = 1
        for size in l1_lists[0].shape:
            norm_array_sizes *= size
        norm_array_sizes *= l1_lists.size

        t_array_size = 1
        for size in t_lists.shape:
            t_array_size *= size

        # print(f"Norm array size : {norm_array_sizes}")
        # print(f"T array sizes : {t_array_size}")

        # SHARING MEMORY BETWEEN PROCESSES
        mp_l1 = mp.Array(ctypes.c_double, norm_array_sizes)
        mp_l2 = mp.Array(ctypes.c_double, norm_array_sizes)
        mp_linf = mp.Array(ctypes.c_double, norm_array_sizes)
        mp_max_divB = mp.Array(ctypes.c_double, t_array_size)
        mp_mean_divB = mp.Array(ctypes.c_double, t_array_size)
        mp_ts = mp.Array(ctypes.c_double, t_array_size)

        l1_arr = np.frombuffer(mp_l1.get_obj())
        l2_arr = np.frombuffer(mp_l2.get_obj())
        linf_arr = np.frombuffer(mp_linf.get_obj())
        max_divB_arr = np.frombuffer(mp_max_divB.get_obj())
        mean_divB_arr = np.frombuffer(mp_mean_divB.get_obj())
        t_arr = np.frombuffer(mp_ts.get_obj())

        l1_arr = l1_arr.reshape(l1_lists.size, l1_lists[0].shape[0], l1_lists[0].shape[1])
        l2_arr = l2_arr.reshape(l1_lists.size, l1_lists[0].shape[0], l1_lists[0].shape[1])
        linf_arr = linf_arr.reshape(l1_lists.size, l1_lists[0].shape[0], l1_lists[0].shape[1])
        max_divB_arr = max_divB_arr.reshape(t_lists.shape)
        mean_divB_arr = mean_divB_arr.reshape(t_lists.shape)
        t_arr = t_arr.reshape(t_lists.shape)

        field_dict_short_list = np.array([dict() for x in range(len(runs))])
        field_dict_list = np.array([dict() for x in range(len(runs))])

        # print(f"L1 arr size : {l1_arr.shape}")
        

        processes = []

        start = time.time()
        id = 0

        for run in runs:
            
            print(f"\n[STARTING] Starting analysis for {run}")

            # Make plot with divB?
            divB = True

            # Get sorted outputs for each run (run denotes set of num. methods)
            if self.runname == "turbulence":
                outputs = sorted(glob.glob(f"{run}/parthenon.prim.?????.{self.athenapksuffix}")) + sorted(glob.glob(f"{run}/Turb.prim.?????.{self.kathenasuffix}"))
                ref_out = sorted(glob.glob(f"{refrun}/Turb.prim.?????.{self.athenapksuffix}")) + sorted(glob.glob(f"{refrun}/Turb.prim.?????.{self.kathenasuffix}"))
            else:
                outputs = sorted(glob.glob(f"{run}/*.{self.athenapksuffix}")) + sorted(glob.glob(f"{run}/*.{self.kathenasuffix}"))
                ref_out = sorted(glob.glob(f"{refrun}/*.{self.athenapksuffix}")) + sorted(glob.glob(f"{refrun}/*.{self.kathenasuffix}"))
           
            # Get num methods used for run and store it for later labelling
            print(f"Sorted athenapk suffix = {sorted(glob.glob(f'{run}/Turb.prim.?????.{self.athenapksuffix}'))}")
            riemann = outputs[0].split("/")[-2].split("_")[-2]
            reconst = outputs[0].split("/")[-2].split("_")[-1]
            code    = outputs[0].split("/")[-2].split("_")[-4]
            integ   = outputs[0].split("/")[-2].split("_")[-3]

            print(f"refrun = {refrun}")

            self.refcode = ref_out[0].split("/")[-2].split("_")[-4]

            print(id, outputs[0], code)

            print(f"ACTUAL RIEMANN : {riemann}")
            riemanns[id] = riemann
            recons[id] = reconst
            codes[id] = code
            integs[id] = integ
            print(f"RIEMANN IN LIST : {riemanns[id]}")
            print(f"CODE APPENDED TO LIST : {codes[id]}, id : {id}")

            if codes[id] == "athenapk":
                field_dict_list[id]       = self.athenapk_field_dict
                field_dict_short_list[id] = self.athenapk_field_dict_short
            elif codes[id] == "kathena":
                field_dict_list[id]       = self.kathena_field_dict
                field_dict_short_list[id] = self.kathena_field_dict_short
            else:
                print("ERROR, unrecognizable code when setting up field dicts...")


            # Prepare each process to be run
            processes.append(mp.Process(target=self.runThread, args=(id, codes, divB, run, refrun, 
                                                                    l1_arr, l2_arr, 
                                                                    linf_arr, max_divB_arr, mean_divB_arr, t_arr, field_dict_list, field_dict_short_list)))
            print(f"L1 lists {id}: {l1_lists[id]}")
            
            id += 1
        
        id = 0
        for p in processes:
            print(f"\n[PARALLEL] Spawning thread {id} for {runs[id]}")
            p.start()
            id += 1

        id = 0
        for p in processes:
            p.join()
            print(f"\n[PARALLEL] Joining thread {id} for {runs[id]}")
            id += 1

        # print(f"L1 lists after : {l1_arr}")
        print(f"Riemanns: {riemanns}\n\nRecons: {recons}")

        if self.saveim == True:
            self.getNorms(l1_arr, l2_arr, linf_arr, t_arr, riemanns, recons, codes, field_list)

            self.getDivBPlot(max_divB_arr, mean_divB_arr, t_arr, recons, riemanns, codes)

        self.writeDataToOutfile(runs, l1_arr, l2_arr, linf_arr, max_divB_arr, mean_divB_arr, t_arr)

        end = time.time()

        print(f"[FINISHED] Finished analysis in {end-start} sec. ({(end-start)/60} min.)")

        return 0




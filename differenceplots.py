import yt
import glob
import numpy as np 
import matplotlib as mpl 
mpl.use("Agg")
import matplotlib.pyplot as plt 
import sys

if __name__ == "__main__":

    runname = sys.argv[1]
    refname = sys.argv[2]

    run_outputs = sorted(glob.glob(f"data/{runname}/*.phdf"))
    ref_outputs = sorted(glob.glob(f"data/{refname}/*.phdf"))

    print(len(run_outputs))
    for i in range(len(run_outputs)):
        print(f"Making plot for dataset {run_outputs[i]}")
        ds     = yt.load(run_outputs[i])
        ref_ds = yt.load(ref_outputs[i])

        ad = ds.all_data()
        ref_ad = ref_ds.all_data()

        num_fields = len(ds.field_list)
        fig, ax = plt.subplots(3, num_fields, figsize=(10,5))

        for (j,field) in enumerate(ds.field_list):

            arr = ad[field]
            ref_arr = ref_ad[field]

            data = arr.to_ndarray().reshape(ds.domain_dimensions)
            ref_data = ref_arr.to_ndarray().reshape(ref_ds.domain_dimensions)

            diff = ref_data - data

            ax[0,j].imshow(ref_data[:,:,0])
            ax[1,j].imshow(data[:,:,0])
            ax[2,j].imshow(diff[:,:,0], cmap="jet")

            ax[0,j].set_title(field[-1])

        ax[0,0].set_ylabel("Reference")
        ax[1,0].set_ylabel("Comparison")
        ax[2,0].set_ylabel("Magnitude Diff.")

        fig.text(0.5, 0.0001, s=f"")

        # fig.tight_layout()

        print(f"Saving to : {runname}_{str(i).zfill(5)}_diff_plot.png")
        plt.savefig(f"{runname}_{str(i).zfill(5)}_diff_plot.png",dpi=276)
        plt.close()

    # USED FOR COMPARING ACROSS SIMULATIONS
    # for (i,sim) in enumerate(sims):

    #     outputs.append(glob.glob(f"{sim}/*.phdf"))


    # for i in range(len(outputs[0])):

    #     ds_list = []

    #     for j in range(len(outputs)):
    #         ds = yt.load(outputs[j][i])

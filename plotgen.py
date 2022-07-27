from cmath import e
import numpy as np 
import matplotlib.pyplot as plt
import h5py
import palettable
import os

# the style_dict key corresponds to the directory name in which the analysis data is stored
sim_dict = {
    'REF_turbulence_10T_sub_kathena_rk4_hlld_4' : {
        'style' : {
            'label' : 'FID-CT_RK4_HLLD_PPM+T',
            'color' : "firebrick",
        },
        'first_dump_id' : 1,  # id (int) of first dump to be read
        'last_dump_id' : 101,  # id (int) of last dump to be read
        'analysis_res' : 256, # resolution that was used for running the analysis
        # *dynamical or turnover* time between dumps
        # assuming equally spaced (in time) data dumps
        'time_between_dumps' : 0.1,
        # the filename pattern should include named variables for 
        # "full_data_dir" (directory), 
        # "dump" (the dump number), and
        # "res", the analysis resolution
        # 'filename_pattern' : '{full_data_dir:s}/{dump:04d}-stats-{res:d}.hdf5',
        'filename_pattern' : '{full_data_dir:s}/flow_analysis_10T_sub_rk4_hlld_4_{dump:05d}.hdf5',
    },

    'turbulence_10T_sub_kathena_rk4_hlle_4' : {
        'style' : {
            'label' : 'CT_RK4_HLLE_PPM+T',
            'color' : palettable.colorbrewer.qualitative.Set1_9.mpl_colors[1],
        },
        'first_dump_id' : 1,  
        'last_dump_id' : 101,  
        'analysis_res' : 256, 
        'time_between_dumps' : 0.1,
        'filename_pattern' : '{full_data_dir:s}/flow_analysis_10T_sub_rk4_hlle_4_{dump:05d}.hdf5',
    },

    'turbulence_10T_sub_kathena_rk4_roe_4' : {
        'style' : {
            'label' : 'CT_RK4_ROE_PPM+T',
            'color' : palettable.colorbrewer.qualitative.Set1_9.mpl_colors[2],
        },
        'first_dump_id' : 1,  
        'last_dump_id' : 101,  
        'analysis_res' : 256, 
        'time_between_dumps' : 0.1,
        'filename_pattern' : '{full_data_dir:s}/flow_analysis_10T_sub_rk4_roe_4_{dump:05d}.hdf5',
    },

    'turbulence_10T_sub_kathena_rk3_hlld_4' : {
        'style' : {
            'label' : 'CT_RK3_HLLD_PPM+T',
            'color' : palettable.colorbrewer.qualitative.Set1_9.mpl_colors[3],
        },
        'first_dump_id' : 1,  
        'last_dump_id' : 101,  
        'analysis_res' : 256, 
        'time_between_dumps' : 0.1,
        'filename_pattern' : '{full_data_dir:s}/flow_analysis_10T_sub_rk3_hlld_4_{dump:05d}.hdf5',
    },

    'turbulence_10T_sub_kathena_rk2_hlld_4' : {
        'style' : {
            'label' : 'CT_RK2_HLLD_PPM+T',
            'color' : palettable.colorbrewer.qualitative.Set1_9.mpl_colors[4],
        },
        'first_dump_id' : 1,  
        'last_dump_id' : 101,  
        'analysis_res' : 256, 
        'time_between_dumps' : 0.1,
        'filename_pattern' : '{full_data_dir:s}/flow_analysis_10T_sub_rk2_hlld_4_{dump:05d}.hdf5',
    },

    'turbulence_10T_sub_kathena_vl2_hlld_4' : {
        'style' : {
            'label' : 'CT_VL2_HLLD_PPM+T',
            'color' : palettable.colorbrewer.qualitative.Set1_9.mpl_colors[5],
        },
        'first_dump_id' : 1,  
        'last_dump_id' : 101,  
        'analysis_res' : 256, 
        'time_between_dumps' : 0.1,
        'filename_pattern' : '{full_data_dir:s}/flow_analysis_10T_sub_vl2_hlld_4_{dump:05d}.hdf5',
    },

    'turbulence_10T_sub_kathena_rk4_hlld_3' : {
        'style' : {
            'label' : 'CT_RK4_HLLD_PPM',
            'color' : palettable.colorbrewer.qualitative.Set1_9.mpl_colors[6],
        },
        'first_dump_id' : 1,  
        'last_dump_id' : 101,  
        'analysis_res' : 256, 
        'time_between_dumps' : 0.1,
        'filename_pattern' : '{full_data_dir:s}/flow_analysis_10T_sub_rk4_hlld_3_{dump:05d}.hdf5',
    },

    'turbulence_10T_sub_kathena_rk4_hlld_2' : {
        'style' : {
            'label' : 'CT_RK4_HLLD_PLM',
            'color' : palettable.colorbrewer.qualitative.Set1_9.mpl_colors[7],
        },
        'first_dump_id' : 1,  
        'last_dump_id' : 101,  
        'analysis_res' : 256, 
        'time_between_dumps' : 0.1,
        'filename_pattern' : '{full_data_dir:s}/flow_analysis_10T_sub_rk4_hlld_2_{dump:05d}.hdf5',
    },

    'turbulence_10T_sub_athenapk_rk3_hlld_ppm' : {
        'style' : {
            'label' : 'DC_RK3_HLLD_PPM',
            'color' : palettable.colorbrewer.qualitative.Set1_9.mpl_colors[8],
        },
        'first_dump_id' : 1,  
        'last_dump_id' : 99,  
        'analysis_res' : 256, 
        'time_between_dumps' : 0.1,
        'filename_pattern' : '{full_data_dir:s}/flow_analysis_10T_athenapk_sub_rk3_hlld_ppm_{dump:05d}.hdf5',
    },
}

grouping_dict = {
    "integrators" : [
        'REF_turbulence_10T_sub_kathena_rk4_hlld_4',
        'turbulence_10T_sub_kathena_rk3_hlld_4',
        'turbulence_10T_sub_kathena_rk2_hlld_4',
        'turbulence_10T_sub_kathena_vl2_hlld_4'
    ],

    "riemanns" : [
        'REF_turbulence_10T_sub_kathena_rk4_hlld_4',
        'turbulence_10T_sub_kathena_rk4_hlle_4',
        'turbulence_10T_sub_kathena_rk4_roe_4'
    ],

    "reconstruction" : [
        'REF_turbulence_10T_sub_kathena_rk4_hlld_4',
        'turbulence_10T_sub_kathena_rk4_hlld_3',
        'turbulence_10T_sub_kathena_rk4_hlld_2'
    ],

    "divergence" : [
        "turbulence_10T_sub_kathena_rk3_hlld_4",
        "turbulence_10T_sub_athenapk_rk3_hlld_ppm"
    ]
}

Alpha = 0.3 # Alpha of the lines not currently being emphasize in panels
Discard = 3.0 # How many turnover times to discard at start of sim for averaging
Low = 7
Up = 40

Comp = {
    "rhoU" : 5./3.,
    "B" : 1.7,
    "u" : 5./3.
}


def fetchData(FlowStats, sim_dict, Ids, print_miss_files=False):

    print_missing_files_warn = print_miss_files
    for Id in Ids:
        
        FlowStats[Id] = {}
            
        AnaRes = sim_dict[Id]['analysis_res']
            
        FlowStats[Id][AnaRes] = {}
        readFiles = 0         

        for Dump in np.arange(sim_dict[Id]['first_dump_id'],
                            sim_dict[Id]['last_dump_id'] + 1):
            
            if Id.split("_")[-4] == "athenapk":
                this_dir = AthenaRootDir + Id
            else:
                this_dir = RootDir + Id

            # print(f"this_dir = {this_dir}")

            this_file = sim_dict[Id]['filename_pattern'].format(
                full_data_dir = this_dir,
                dump = Dump,
                res = AnaRes
            )
            if not os.path.isfile(this_file):
                if print_missing_files_warn:
                    print("missing " + this_file)
                continue

            FlowStats[Id][AnaRes][Dump] = {}

            try:
                flowquant = h5py.File(this_file,"r")

                FlowStats[Id][AnaRes][Dump] = flowquant                

                readFiles += 1

            except:                    
                del FlowStats[Id][AnaRes][Dump]
                print("bad file " + this_file)


        print("Got %d flowQuants for %s at res %d" % (readFiles,Id,AnaRes))

def getTemporalData(Id, AnaRes, QuanName, FlowStats):

    tmpY = []
    tmpX = []
    for Dump in sorted(FlowStats[Id][AnaRes].keys()):
        if Dump * sim_dict[Id]["time_between_dumps"] < Discard:
            continue
        tmpX.append(sim_dict[Id]['time_between_dumps']*Dump)
        tmpY.append(FlowStats[Id][AnaRes][Dump][QuanName])

    return tmpX, tmpY

def plotTimeAveragedQuant(QuanStruct, FlowStats, sim_dict, fname):

    QuanName = QuanStruct[0]
    QuanLabel = QuanStruct[1]

    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6,10))

    # Loop over the sims from changing integrator and plot
    for i, sim in enumerate(grouping_dict["integrators"]):
        time, quant = getTemporalData(sim, sim_dict[sim]["analysis_res"], QuanName, FlowStats)
        ax[0].plot(time, quant, **sim_dict[sim]["style"])

    # Loop over the sims from changing riemann and plot
    for i, sim in enumerate(grouping_dict["riemanns"]):
        time, quant = getTemporalData(sim, sim_dict[sim]["analysis_res"], QuanName, FlowStats)
        ax[1].plot(time, quant, **sim_dict[sim]["style"])

    # Loop over the sims from changing reconstruction and plot
    for i, sim in enumerate(grouping_dict["reconstruction"]):
        time, quant = getTemporalData(sim, sim_dict[sim]["analysis_res"], QuanName, FlowStats)
        ax[2].plot(time, quant, **sim_dict[sim]["style"])

    # Loop over everything, add labels, and grayscale lines in back
    ax[2].set_xlabel("Time")
    for i in range(3):
        if i == 0:
            type = "by Integrator"
        elif i == 1:
            type = "by Riemann"
        else:
            type = "by Reconst."

        ax[i].set_ylabel(QuanLabel + "\n" + type)
        ax[i].legend()
        ax[i].grid()

        for Id in list(sim_dict.keys()):
            time, quant = getTemporalData(Id, sim_dict[Id]["analysis_res"], QuanName, FlowStats)
            ax[i].plot(time, quant, color="black", alpha=Alpha, zorder=1)

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.00,hspace=0.01)

    plt.savefig(fname)
    plt.close()

    # Comparing Codes
    fig, ax = plt.subplots(1, 1)

    for i, sim in enumerate(grouping_dict["divergence"]):
        time, quant = getTemporalData(sim, sim_dict[sim]["analysis_res"], QuanName, FlowStats)
        ax.plot(time, quant, **sim_dict[sim]["style"])

    ax.legend()
    ax.grid()

    ax.set_ylabel(QuanLabel + "\n" + "DC vs CT")
    ax.set_xlabel("Time")

    fig.tight_layout()

    plt.savefig("comparing_codes_" + QuanStruct[2] + ".png")
    plt.close()

    return

def getTimeAveragedSpectra(FlowStats, Id, AnaRes, Quan, Type):
    tmp = []    
    for Dump in FlowStats[Id][AnaRes].keys():
        if Dump * sim_dict[Id]['time_between_dumps'] < Discard:
            continue
        tmp.append(FlowStats[Id][AnaRes][Dump][Quan + '/PowSpec/'+Type][1])
    return FlowStats[Id][AnaRes][Dump][Quan + '/PowSpec/'+Type][0],  np.mean(tmp,axis=0), np.std(tmp,axis=0)
    
def plotAveragedSpectra(FlowStats, sim_dict, QuanStruct, fname, show_std=False):

    Quan = QuanStruct[0]
    Type = QuanStruct[1]
    QuanLabel = QuanStruct[2]

    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6,10))

    # Loop over the sims from changing integrator and plot
    for i, sim in enumerate(grouping_dict["integrators"]):
        AnaRes = sim_dict[sim]["analysis_res"]
        X, Y, Yerr = getTimeAveragedSpectra(FlowStats, sim, AnaRes, Quan, Type) # Type can be "Full", "Sol", or "Dil"
        
        mask = np.logical_and(X < AnaRes/2, X!=0.0)
        
        ax[0].plot(X[mask], Y[mask], **sim_dict[sim]["style"])

        if show_std:
            ax[0].fill_between(X[mask], Y[mask]-Yerr[mask], Y[mask]+Yerr[mask], alpha=0.3, color=sim_dict[sim]["style"]["color"])

    # Loop over the sims from changing riemann and plot
    for i, sim in enumerate(grouping_dict["riemanns"]):
        AnaRes = sim_dict[sim]["analysis_res"]
        X, Y, Yerr = getTimeAveragedSpectra(FlowStats, sim, AnaRes, Quan, Type) # Type can be "Full", "Sol", or "Dil"
        
        mask = np.logical_and(X < AnaRes/2, X!=0.0)
        
        ax[1].plot(X[mask], Y[mask], **sim_dict[sim]["style"])

        if show_std:
            ax[1].fill_between(X[mask], Y[mask]-Yerr[mask], Y[mask]+Yerr[mask], alpha=0.3, color=sim_dict[sim]["style"]["color"])

    # Loop over the sims from changing reconstruction and plot
    for i, sim in enumerate(grouping_dict["reconstruction"]):
        AnaRes = sim_dict[sim]["analysis_res"]
        X, Y, Yerr = getTimeAveragedSpectra(FlowStats, sim, AnaRes, Quan, Type) # Type can be "Full", "Sol", or "Dil"
        
        mask = np.logical_and(X < AnaRes/2, X!=0.0)
        
        ax[2].plot(X[mask], Y[mask], **sim_dict[sim]["style"])

        if show_std:
            ax[2].fill_between(X[mask], Y[mask]-Yerr[mask], Y[mask]+Yerr[mask], alpha=0.3, color=sim_dict[sim]["style"]["color"])

    ax[2].set_xlabel("Wavenumber k")
    for i in range(3):
        if i == 0:
            type = "by Integrator"
        elif i == 1:
            type = "by Riemann"
        else:
            type = "by Reconst."

        ax[i].set_ylabel(QuanLabel + "\n" + type)
        ax[i].legend()
        ax[i].grid()

        for Id in list(sim_dict.keys()):
            X, Y, Yerr = getTimeAveragedSpectra(FlowStats, Id, AnaRes, Quan, Type)

            mask = np.logical_and(X < AnaRes/2, X!=0.0)

            ax[i].plot(X[mask], Y[mask], color="black", alpha=Alpha, zorder=1)

    for i in range(3):
        if Quan != "B":
            ax[i].axvspan(Low,Up, facecolor='black', alpha=0.05)
        ax[i].set_xlim(1.,AnaRes/2)
        ax[i].set_xscale("log")
        ax[i].set_yscale("log")

    ax[0].tick_params(axis='x', direction='in')
    ax[0].tick_params(axis='x',tickdir='in')
    ax[1].tick_params(axis='x', direction='in')
    ax[1].tick_params(axis='x',tickdir='in')

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.00,hspace=0.01)

    plt.savefig(fname)
    plt.close()

    # Comparing codes
    fig, ax = plt.subplots(1, 1)

    for i, sim in enumerate(grouping_dict["divergence"]):
        AnaRes = sim_dict[sim]["analysis_res"]
        X, Y, Yerr = getTimeAveragedSpectra(FlowStats, sim, AnaRes, Quan, Type) # Type can be "Full", "Sol", or "Dil"
        
        mask = np.logical_and(X < AnaRes/2, X!=0.0)
        
        ax.plot(X[mask], Y[mask], **sim_dict[sim]["style"])

        if show_std:
            ax.fill_between(X[mask], Y[mask]-Yerr[mask], Y[mask]+Yerr[mask], alpha=0.3, color=sim_dict[sim]["style"]["color"])

    ax.legend()
    ax.grid()
    ax.set_xlim(1.,AnaRes/2)
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_ylabel(QuanLabel + "\n" + "DC vs CT")
    fig.tight_layout()

    plt.savefig("comparing_codes_spectra_" + Quan + ".png")
    plt.close()

    return

def plotCompensatedAveragedSpectra(FlowStats, sim_dict, QuanStruct, Comp, fname, show_std=False):

    Quan = QuanStruct[0]
    Type = QuanStruct[1]
    QuanLabel = QuanStruct[2]

    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6,10))

    # Loop over the sims from changing integrator and plot
    for i, sim in enumerate(grouping_dict["integrators"]):
        AnaRes = sim_dict[sim]["analysis_res"]
        X, Y, Yerr = getTimeAveragedSpectra(FlowStats, sim, AnaRes, Quan, Type) # Type can be "Full", "Sol", or "Dil"
        
        mask = np.logical_and(X < AnaRes/2, X!=0.0)
        
        ax[0].plot(X[mask], X[mask]**Comp[Quan]*Y[mask], **sim_dict[sim]["style"])

        if show_std:
            ax[0].fill_between(X[mask], X[mask]**Comp[Quan] *(Y[mask]-Yerr[mask]), 
                               X[mask]**Comp[Quan]*(Y[mask]+Yerr[mask]), alpha=0.3, 
                               color=sim_dict[sim]["style"]["color"])

    # Loop over the sims from changing riemann and plot
    for i, sim in enumerate(grouping_dict["riemanns"]):
        AnaRes = sim_dict[sim]["analysis_res"]
        X, Y, Yerr = getTimeAveragedSpectra(FlowStats, sim, AnaRes, Quan, Type) # Type can be "Full", "Sol", or "Dil"
        
        mask = np.logical_and(X < AnaRes/2, X!=0.0)
        
        ax[1].plot(X[mask], X[mask]**Comp[Quan]*Y[mask], **sim_dict[sim]["style"])

        if show_std:
            ax[1].fill_between(X[mask], X[mask]**Comp[Quan] *(Y[mask]-Yerr[mask]), 
                               X[mask]**Comp[Quan] *(Y[mask]+Yerr[mask]), alpha=0.3, 
                               color=sim_dict[sim]["style"]["color"])

    # Loop over the sims from changing reconstruction and plot
    for i, sim in enumerate(grouping_dict["reconstruction"]):
        AnaRes = sim_dict[sim]["analysis_res"]
        X, Y, Yerr = getTimeAveragedSpectra(FlowStats, sim, AnaRes, Quan, Type) # Type can be "Full", "Sol", or "Dil"
        
        mask = np.logical_and(X < AnaRes/2, X!=0.0)
        
        ax[2].plot(X[mask], X[mask]**Comp[Quan]*Y[mask], **sim_dict[sim]["style"])

        if show_std:
            ax[2].fill_between(X[mask], X[mask]**Comp[Quan] *(Y[mask]-Yerr[mask]), 
                               X[mask]**Comp[Quan] *(Y[mask]+Yerr[mask]), alpha=0.3, 
                               color=sim_dict[sim]["style"]["color"])

    ax[2].set_xlabel("Wavenumber k")
    for i in range(3):
            
        if i == 0:
            type = "by Integrator"
        elif i == 1:
            type = "by Riemann"
        else:
            type = "by Reconst."

        ax[i].set_ylabel("Comp. " + QuanLabel + "\n" + type)
        ax[i].legend()
        ax[i].grid()

        for Id in list(sim_dict.keys()):
            X, Y, Yerr = getTimeAveragedSpectra(FlowStats, Id, AnaRes, Quan, Type)

            mask = np.logical_and(X < AnaRes/2, X!=0.0)

            ax[i].plot(X[mask], X[mask]**Comp[Quan]*Y[mask], color="black", alpha=Alpha, zorder=1)

    for i in range(3):
        if Quan != "B":
            ax[i].axvspan(Low,Up, facecolor='black', alpha=0.05)
        ax[i].set_xlim(1.,AnaRes/2)
        ax[i].set_xscale("log")
        ax[i].set_yscale("log")

    ax[0].tick_params(axis='x', direction='in')
    ax[0].tick_params(axis='x',tickdir='in')
    ax[1].tick_params(axis='x', direction='in')
    ax[1].tick_params(axis='x',tickdir='in')

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.00,hspace=0.01)

    plt.savefig(fname)
    plt.close()

    # Comparing codes
    fig, ax = plt.subplots(1, 1)

    for i, sim in enumerate(grouping_dict["divergence"]):
        AnaRes = sim_dict[sim]["analysis_res"]
        X, Y, Yerr = getTimeAveragedSpectra(FlowStats, sim, AnaRes, Quan, Type) # Type can be "Full", "Sol", or "Dil"
        
        mask = np.logical_and(X < AnaRes/2, X!=0.0)
        
        ax.plot(X[mask], X[mask]**Comp[Quan]*Y[mask], **sim_dict[sim]["style"])

        if show_std:
            ax.fill_between(X[mask], X[mask]**Comp[Quan] *(Y[mask]-Yerr[mask]), 
                               X[mask]**Comp[Quan] *(Y[mask]+Yerr[mask]), alpha=0.3, 
                               color=sim_dict[sim]["style"]["color"])
    ax.legend()
    ax.grid()
    ax.set_xlim(1.,AnaRes/2)
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_ylabel("Comp. " + QuanLabel + "\n" + "DC vs CT")
    fig.tight_layout()

    plt.savefig("comparing_codes_comp_spectra_" + Quan + ".png")
    plt.close()

    return

if __name__ == "__main__":

    # root directory of the simulation
    # This currently assumes that all subdirectories (= style_dict keys) are in this single root directory
    RootDir = '/mnt/scratch/fushstep/turb-sims/'
    AthenaRootDir = "/mnt/home/fushstep/athenapk/turb-sims/"
    PlotPrefix = "plot"

    Ids = list(sim_dict.keys())

    FlowStats = {}

    # Getting the data for each of the simulations
    fetchData(FlowStats, sim_dict, Ids)

    quantities_to_plot = [
        ('u' + '/moments/' + 'rms',r'RMS $\mathrm{Ms}$',"rmsU"),
        ('KinEnSpecific' + '/moments/' + 'mean',r'Mean $\mathrm{KE_s}$', "meanKEspec"),
        ('KinEnDensity' + '/moments/' + 'mean',r'Mean $\mathrm{KE}$', "meanKEdens"),
        ('AlfvenicMach' + '/moments/' + 'rms',r'RMS $\mathrm{M_a}$', "Ma"),                                  
        ('rho-B/corr',r'$\mathrm{Corr}[\rho,B]$', "rho-Bcorr"),
        ('MagEnDensity' + '/moments/' + 'mean', r'Mean $E_B$', "magEn")
    ]

    for QuanStruct in quantities_to_plot:
        fname = f"{PlotPrefix}_{QuanStruct[-1]}.png"
        print(f"Making temporal evolution plot ({fname}) for {QuanStruct[0]}")
        plotTimeAveragedQuant(QuanStruct, FlowStats, sim_dict, fname)


    sprectra_to_plot = [
        ("rhoU", "Full", "Kin. Energy Spectrum"),
        ("B", "Full", "Magnetic Energy Spectrum")
    ]

    for QuanStruct in sprectra_to_plot:
        fname = f"{PlotPrefix}_spectra_nostd_{QuanStruct[0]}.png"
        print(f"Making {QuanStruct[1]} temporal averaged spectra plot ({fname}) for {QuanStruct[0]}")
        plotAveragedSpectra(FlowStats, sim_dict, QuanStruct, fname, show_std=False)

        fname = f"{PlotPrefix}_comp_spectra_nostd_{QuanStruct[0]}.png"
        print(f"Making {QuanStruct[1]} compensated temporal averaged spectra plot ({fname}) for {QuanStruct[0]}")
        plotCompensatedAveragedSpectra(FlowStats, sim_dict, QuanStruct, Comp, fname, show_std=False)

    for QuanStruct in sprectra_to_plot:
        fname = f"{PlotPrefix}_spectra_std_{QuanStruct[0]}.png"
        print(f"Making {QuanStruct[1]} temporal averaged spectra plot ({fname}) for {QuanStruct[0]}")
        plotAveragedSpectra(FlowStats, sim_dict, QuanStruct, fname, show_std=True)

        fname = f"{PlotPrefix}_comp_spectra_std_{QuanStruct[0]}.png"
        print(f"Making {QuanStruct[1]} compensated temporal averaged spectra plot ({fname}) for {QuanStruct[0]}")
        plotCompensatedAveragedSpectra(FlowStats, sim_dict, QuanStruct, Comp, fname, show_std=True)

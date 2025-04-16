import os
from os import path as p
import matplotlib
matplotlib.use('Agg')  # Set the Agg backend before importing pyplot
import matplotlib.pyplot as plt
import numpy as np
from subprocess import call, PIPE
import imageio
#####################################################################
def extract_number(filename):
    # Extract the number from the filename
    return int(filename.split('_')[-1].split('.')[0])
    
#####################################################################

def make_gif(inDir, outGif):

    pngFiles = sorted(
    [file for file in os.listdir(inDir) 
     if "fitting_shuffle" in file and file.endswith(".png")],
    key=extract_number)    

    images = []
    for file in pngFiles:
        images.append(imageio.v3.imread(p.join(inDir, file)))

    imageio.mimwrite(outGif, images, 'GIF', fps=2, loop = 0)



#####################################################################
def set_rc_params():
    darkGrey: str = '#1a1a1a'

    yellow: str = '#FFFF00'

    ## set text and axes colors
    plt.rcParams['text.color'] = yellow
    plt.rcParams['axes.labelcolor'] = yellow
    plt.rcParams['xtick.color'] = yellow
    plt.rcParams['ytick.color'] = yellow
    plt.rcParams['axes.edgecolor'] = yellow
    ## set the background color
    plt.rcParams['axes.facecolor'] = darkGrey
    plt.rcParams['figure.facecolor'] = darkGrey
    ## sort out font sizes
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.titlesize'] = 22
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['legend.fontsize'] = 18



#####################################################################
def plot_qmmm_energies(qmTotalEnergy,
                        qmTorsionEnergy,
                        mmTotalEnergy, 
                        mmTorsionEnergy,
                        cosineComponents,
                        outDir,
                        shuffleIndex):
    

    set_rc_params()
    ## init some colors to be used
    white :str = '#FFFFFF'
    magenta : str = '#FF00FF'
    brightCyan: str = '#00FFFF'
    ## construct angles
    angles = np.linspace(0, 360, len(qmTotalEnergy))

    fig, axis = plt.subplots( nrows=1, ncols=2, figsize=(12, 8))  # create figure
    fig.subplots_adjust(bottom=0.25)  # Adjust the bottom margin

    ## plot QM and MM total energies GLOW
    for n in range(1, 7):
        axis[0].plot(angles, qmTotalEnergy, color=magenta, linewidth=1+n, alpha=0.1)
        axis[0].plot(angles, mmTotalEnergy, color=brightCyan, linewidth=1+n, alpha=0.1)
    ## plot main traces
    axis[0].plot(angles, qmTotalEnergy, label='QM', linewidth=2, color=magenta)
    axis[0].plot(angles, mmTotalEnergy, label='MM', linewidth=2, color=brightCyan)
    ## LEGEND ##
    axis[0].legend(['QM', 'MM'],
                    loc="upper center",
                    bbox_to_anchor=(0.5, -0.12),
                    ncol=2,
                    handlelength=0,
                    handletextpad=0,
                    labelcolor=[magenta, brightCyan])
    
    axis[0].set_xlabel('Torsion Angle')
    axis[0].set_ylabel('Energy (Kcal / mol)')
    axis[0].set_title("Total Energies")
    axis[0].set_xlim(0, 360)
    axis[0].set_xticks(np.arange(0, 361, 60))

    ## plot QM and MM torsion energies GLOW
    for n in range(1, 7):
        axis[1].plot(angles, qmTorsionEnergy, color=magenta, linewidth=1+n, alpha=0.1)
        axis[1].plot(angles, mmTorsionEnergy, color=brightCyan, linewidth=1+n, alpha=0.1)
    ## plot main traces
    axis[1].plot(angles, qmTorsionEnergy, label='QM', linewidth=2, color=magenta)
    axis[1].plot(angles, mmTorsionEnergy, label='MM', linewidth=2, color=brightCyan)

    ## plot cosine components
    for freq, cosine_component in cosineComponents.items():
        axis[1].plot(angles, cosine_component, linestyle='--',
                     label=f'Cosine Component {freq:.0f}', color=white, alpha=0.5)
    # Dummy plot for legend
    axis[1].plot([], [], linestyle='--', color=white, alpha=0.5,
                 label='Parameters')
    ## LEGEND ##
    axis[1].legend(['QM', 'MM', 'Parameters'],
                    loc="upper center",
                    bbox_to_anchor=(0.5, -0.12),
                    ncol=3,
                    handlelength=0,
                    handletextpad=0,
                    labelcolor=[magenta, brightCyan, white])

    axis[1].set_xlabel('Torison Angle')
    axis[1].set_ylabel('Energy (Kcal / mol)')
    axis[1].set_title("Torsion Energies")
    axis[1].set_xlim(0, 360)
    axis[1].set_xticks(np.arange(0, 361, 60))

    ## add a shuffle index label
    fig.text(0.1, 0.9, f'{shuffleIndex + 1}', fontsize=16,
             color='yellow', bbox=dict(facecolor='none', 
             edgecolor='yellow', boxstyle='round,pad=0.5'))


    fig.savefig(p.join(outDir, f"fitting_shuffle_{shuffleIndex+1}.png"))
    plt.close(fig)
#####################################################################


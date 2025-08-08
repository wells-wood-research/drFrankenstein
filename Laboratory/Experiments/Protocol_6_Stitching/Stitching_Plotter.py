import os
from os import path as p
import matplotlib
matplotlib.use('Agg')  # Set the Agg backend before importing pyplot
import matplotlib.pyplot as plt
import numpy as np
from subprocess import call, PIPE
import os
import pandas as pd
from PIL import Image


## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass

def plot_run_mean_average_error(outDir, maeTorsionDf: pd.DataFrame, maeTotalDf: pd.DataFrame):
    set_rc_params()
    ## init some colors to be used
    white: str = '#FFFFFF'
    brightGreen: str = '#00FF00'
    brightOrange: str = '#FFA500'
    fig, axis = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))  # create figure

    ## plot average MAE with GLOW
    for n in range(1, 7):
        axis.plot(maeTorsionDf["All_Torsions"], color=brightGreen, linewidth=1+n*1.25, alpha=0.1)
    ## plot main traces
    axis.plot(maeTorsionDf["All_Torsions"], label='Torsion', linewidth=2, color=brightGreen)
    ## plot subtle traces for individual torsions
    individualTorsionDf = maeTorsionDf.drop(["All_Torsions"], axis=1)
    for torsionTag in individualTorsionDf.columns:
        axis.plot(individualTorsionDf[torsionTag], linewidth=1, color=brightGreen, alpha=0.3)

    ## plot average MAE with GLOW
    for n in range(1, 7):
        axis.plot(maeTotalDf["All_Torsions"], color=brightOrange, linewidth=1+n*1.25, alpha=0.1)
    ## plot main traces
    axis.plot(maeTotalDf["All_Torsions"], label='Total', linewidth=2, color=brightOrange)
    ## plot subtle traces for individual torsions
    individualTorsionDf = maeTotalDf.drop(["All_Torsions"], axis=1)
    for torsionTag in individualTorsionDf.columns:
        axis.plot(individualTorsionDf[torsionTag], linewidth=1, color=brightOrange, alpha=0.3)

    # Add labels and legend
    axis.set_xlabel('Number of Fitting Shuffles')
    axis.set_ylabel('Mean Average Error (Kcal / mol)')
    axis.set_title('Mean Average Error vs Number of Fitting Shuffles')

    # Get handles and labels for the legend
    handles, labels = axis.get_legend_handles_labels()
    # Filter to include only "All_Torsions" entries
    filtered_handles_labels = [(h, l) for h, l in zip(handles, labels) if l in ['Torsion', 'Total']]
    if filtered_handles_labels:
        handles, labels = zip(*filtered_handles_labels)
        axis.legend(handles, labels)
    else:
        axis.legend()

    plt.savefig(p.join(outDir, "run_mean_average_error.png"), bbox_inches='tight')
    plt.close()

def plot_mean_average_error(torsionFittingDir: DirectoryPath, maeCsv: FilePath, torsionTag: str):

    ## get data from csv
    maeDf = pd.read_csv(maeCsv)
    torsionMaeDf = maeDf[maeDf["torsion_tag"] == torsionTag]
    

    set_rc_params()
    ## init some colors to be used
    brightGreen: str = '#00FF00'
    brightOrange: str = '#FFA500'

    fig, axis = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))  # create figure

    ## plot torsion MAE with GLOW
    for n in range(1, 7):
        axis.plot(torsionMaeDf["mae_torsion"], color=brightGreen, linewidth=1+n*1.25, alpha=0.1)
    ## plot main traces
    axis.plot(torsionMaeDf["mae_torsion"], label='Torsion', linewidth=2, color=brightGreen)

    ## plot torsion MAE with GLOW
    for n in range(1, 7):
        axis.plot(torsionMaeDf["mae_total"], color=brightOrange, linewidth=1+n*1.25, alpha=0.1)
    ## plot main traces
    axis.plot(torsionMaeDf["mae_total"], label='Total', linewidth=2, color=brightOrange)

    # Add labels and legend
    axis.set_xlabel('Number of Fitting Shuffles')
    axis.set_ylabel('Mean Average Error (Kcal / mol)')
    axis.set_title('Mean Average Error vs Number of Fitting Shuffles')
    axis.legend()

    plt.savefig(p.join(torsionFittingDir, "mean_average_error.png"), bbox_inches='tight')
    plt.close()
#####################################################################
def extract_number(filename):
    # Extract the number from the filename
    return int(filename.split('_')[-1].split('.')[0])
    
#####################################################################

def make_gif(inDir, outGif):

    """
    Creates a GIF using Pillow directly in a memory-efficient way
    by using a generator expression.
    """
    pngFiles = sorted(
        [os.path.join(inDir, f) for f in os.listdir(inDir)
         if "fitting_shuffle" in f and f.endswith(".png")],
        key=extract_number
    )

    if not pngFiles:
        return
        
    
    try:
        # Open the first image
        img_first = Image.open(pngFiles[0])

        # Create a generator for the remaining images
        # This is the key to low memory usage!
        # Note the parentheses () instead of square brackets [].
        img_generator = (Image.open(f) for f in pngFiles[1:])

        # Save the first image, providing the generator for the rest
        img_first.save(
            fp=outGif,
            format='GIF',
            save_all=True,
            append_images=img_generator,  # Pass the generator here
            duration=500,                 # 500ms per frame = 2 FPS
            loop=0                        # Loop forever
        )

    except Exception as e:
        raise(e)



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
    for freq, cosineComponent in cosineComponents.items():
        axis[1].plot(angles, cosineComponent, linestyle='--',
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


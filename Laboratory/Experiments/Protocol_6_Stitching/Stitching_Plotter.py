import os
from os import path as p
import matplotlib
matplotlib.use('Agg')  # Set the Agg backend before importing pyplot
import matplotlib.pyplot as plt
import numpy as np
import re
from subprocess import call, PIPE
import pandas as pd
from PIL import Image
import tempfile
from itertools import islice
from matplotlib.ticker import MaxNLocator


## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass

def plot_run_mean_average_error(outDir, maeCsv):

    maeDf = pd.read_csv(maeCsv)

    set_rc_params()
    ## init some colors to be used
    white: str = '#FFFFFF'
    brightGreen: str = '#00FF00'
    brightOrange: str = '#FFA500'
    fig, axis = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))  # create figure


    allTorsionsDf = maeDf[maeDf["torsion_tag"] == "All_Torsions"]
    individualTorsionDf = maeDf[maeDf["torsion_tag"] != "All_Torsions"]
    ## ==================== TORSIONS ================================== ##
    ## plot average MAE with GLOW
    for n in range(1, 7):
        axis.plot(allTorsionsDf["shuffle"], allTorsionsDf["mae_torsion"],color=brightGreen, linewidth=1+n*1.25, alpha=0.1)
    ## plot main traces
    axis.plot(allTorsionsDf["shuffle"], allTorsionsDf["mae_torsion"], label='Torsion', linewidth=2, color=brightGreen)

    ## plot subtle traces for individual torsions
    for torsionTag, df in individualTorsionDf.groupby("torsion_tag"):
        axis.plot(df["shuffle"], df["mae_torsion"], linewidth=1, color=brightGreen, alpha=0.4)

    ## ==================== TOTAL ================================== ##
    ## plot average MAE with GLOW
    for n in range(1, 7):
        axis.plot( allTorsionsDf["shuffle"], allTorsionsDf["mae_total"], color=brightOrange, linewidth=1+n*1.25, alpha=0.1)
    ## plot main traces
    axis.plot(allTorsionsDf["shuffle"], allTorsionsDf["mae_total"], label='Total', linewidth=2, color=brightOrange)
    ## plot subtle traces for individual torsions
    for torsionTag, df in individualTorsionDf.groupby("torsion_tag"):
        axis.plot( df["shuffle"], df["mae_total"], linewidth=1, color=brightOrange, alpha=0.4)

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
    # Get unique shuffle values and ensure they are integers
    shuffle_values = allTorsionsDf["shuffle"].unique().astype(int)

    # Set x-ticks, ensuring only values >= 0 and <= shuffle_values[-1]
    axis.set_xticks(shuffle_values)  # Set all integer ticks initially
    axis.xaxis.set_major_locator(MaxNLocator(integer=True, nbins='auto'))  # Auto-adjust ticks
    # Get current ticks and filter to ensure they are >= 0 and <= shuffle_values[-1]
    current_ticks = axis.get_xticks()
    valid_ticks = current_ticks[(current_ticks >= 0) & (current_ticks <= shuffle_values[-1])]
    # Always include the last value if not already present
    if shuffle_values[-1] not in valid_ticks:
        valid_ticks = np.append(valid_ticks, shuffle_values[-1])
    # Sort ticks to maintain order and set them
    axis.set_xticks(np.sort(valid_ticks))

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
def _extract_number(filename):
    # Extract the number from the filename
    return int(filename.split('_')[-1].split('.')[0])


def _extract_fit_frame_key(filename):
    match = re.search(r"fitting_shuffle_(\d+)_nCosines_(\d+)\.png$", os.path.basename(filename))
    if match:
        return int(match.group(2)), int(match.group(1))
    return (0, _extract_number(filename))
    
#####################################################################

def make_gif(inDir: DirectoryPath, outGif: FilePath, batchSize: int = 50, duration: int = 100):
    """
    Creates a GIF from a generator of PNG file paths using a memory-safe,
    pure Python approach.

    Args:
        inDir (str): The directory containing the PNG files.
        outGif (str): The path to save the final GIF file.
        batchSize (int): The number of images to process in each batch.
        duration (int): The duration for each frame in the GIF, in milliseconds.
    """
    pngGenerator = (
        item for item in sorted(
            (os.path.join(inDir, f) for f in os.listdir(inDir) if f.endswith(".png") and f.startswith("fitting_shuffle")),
            key=_extract_fit_frame_key
        )
    )
    frames = [Image.open(png) for png in pngGenerator]

    if not frames:
        print("No PNG files were found to create a GIF.")
        return

    # Save the first frame and append the rest
    frames[0].save(
        outGif,
        save_all=True,
        append_images=frames[1:],  # Append all other frames
        optimize=True,             # Optimize the GIF palette
        duration=duration,
        loop=0
    )
    
    # Clean up the opened image objects
    for frame in frames:
        frame.close()


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
                        mmFittedTorsionEnergy,
                        cosineComponents,
                        torsionMetrics,
                        totalMetrics,
                        maeTorsion,
                        maeTotal,
                        nCosines,
                        converged,
                        outDir,
                        shuffleIndex):
    

    set_rc_params()
    ## init some colors to be used
    white :str = '#FFFFFF'
    magenta : str = '#FF00FF'
    brightCyan: str = '#00FFFF'
    shuffleLabel = f"fitting_shuffle_{shuffleIndex}_nCosines_{nCosines}"
    if not (len(qmTotalEnergy) == len(qmTorsionEnergy) == len(mmTotalEnergy) == len(mmFittedTorsionEnergy)):
        raise ValueError(
            f"{shuffleLabel}: QM/MM arrays must have matching lengths "
            f"(qm_total={len(qmTotalEnergy)}, qm_torsion={len(qmTorsionEnergy)}, "
            f"mm_total={len(mmTotalEnergy)}, mm_torsion={len(mmFittedTorsionEnergy)})."
        )
    ## construct angles on the 10-degree fitting grid used by the scan/protocol
    angleStepDegrees = 10
    angles = np.arange(len(qmTotalEnergy), dtype=float) * angleStepDegrees

    fig, axis = plt.subplots( nrows=1, ncols=2, figsize=(12, 8))  # create figure

    ## plot QM and MM total energies GLOW
    for n in range(1, 7):
        axis[0].plot(angles, qmTotalEnergy, color=magenta, linewidth=1+n, alpha=0.1)
        axis[0].plot(angles, mmTotalEnergy, color=brightCyan, linewidth=1+n, alpha=0.1)
    ## plot main traces
    axis[0].plot(angles, qmTotalEnergy, label='QM target', linewidth=2, color=magenta)
    axis[0].plot(angles, mmTotalEnergy, label='MM fit', linewidth=2, color=brightCyan)
    axis[0].legend(loc="upper right")
    axis[0].text(0.02, 0.98, _format_metrics_box("Total", totalMetrics, maeTotal), transform=axis[0].transAxes,
                 color="red", fontsize=11, va="top", family="monospace",
                 bbox=dict(facecolor="none", edgecolor="red", boxstyle="round,pad=0.3"))

    axis[0].set_xlabel('Torsion Angle')
    axis[0].set_ylabel('Energy (Kcal / mol)')
    axis[0].set_title("Total Energies")
    axis[0].set_xlim(0, 360)
    axis[0].set_xticks(np.arange(0, 361, 60))

    ## plot QM and MM torsion energies GLOW
    for n in range(1, 7):
        axis[1].plot(angles, qmTorsionEnergy, color=magenta, linewidth=1+n, alpha=0.1)
        axis[1].plot(angles, mmFittedTorsionEnergy, color=brightCyan, linewidth=1+n, alpha=0.1)
    ## plot main traces
    axis[1].plot(angles, qmTorsionEnergy, label='QM target', linewidth=2, color=magenta)
    axis[1].plot(angles, mmFittedTorsionEnergy, label='MM fit', linewidth=2, color=brightCyan)
    axis[1].text(0.02, 0.98, _format_metrics_box("Torsion", torsionMetrics, maeTorsion), transform=axis[1].transAxes,
                 color="red", fontsize=11, va="top", family="monospace",
                 bbox=dict(facecolor="none", edgecolor="red", boxstyle="round,pad=0.3"))
    if converged:
        axis[1].text(
            0.72,
            0.08,
            "CONVERGED",
            transform=axis[1].transAxes,
            color="lime",
            fontsize=18,
            fontweight="bold",
            bbox=dict(facecolor="none", edgecolor="lime", boxstyle="round,pad=0.3"),
        )

    ## plot cosine components
    if isinstance(cosineComponents, dict):
        componentItems = cosineComponents.items()
    else:
        componentItems = cosineComponents
    for freq, cosineComponent in componentItems:
        axis[1].plot(angles, cosineComponent, linestyle='--',
                     color=white, alpha=0.5)
    # Dummy plot for legend
    axis[1].plot([], [], linestyle='--', color=white, alpha=0.5,
                 label='Parameters')
    axis[1].legend(loc="upper right")

    axis[1].set_xlabel('Torsion Angle')
    axis[1].set_ylabel('Energy (Kcal / mol)')
    axis[1].set_title("Torsion Energies")
    axis[1].set_xlim(0, 360)
    axis[1].set_xticks(np.arange(0, 361, 60))

    fig.savefig(p.join(outDir, f"{shuffleLabel}.png"))
    plt.close(fig)


def _format_metrics_box(label, metrics, mae):
    return (
        f"{label}\n"
        f"Score {metrics['composite_score']:.3f}\n"
        f"Loc   {metrics['location_score']:.3f}\n"
        f"Amp   {metrics['amplitude_score']:.3f}\n"
        f"Count {metrics['stationary_count_score']:.3f}\n"
        f"nMAE  {metrics['normalized_mae_score']:.3f}\n"
        f"MAE   {mae:.3f}"
    )
#####################################################################

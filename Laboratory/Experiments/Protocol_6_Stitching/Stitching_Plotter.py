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

import matplotlib.gridspec as gridspec

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
    plt.rcParams.update({'font.family': 'monospace'}) 


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
                        shuffleIndex,
                        tol):
    
    set_rc_params()
    white :str = '#FFFFFF'
    magenta : str = '#FF00FF'
    brightCyan: str = '#00FFFF'
    shuffleLabel = f"fitting_shuffle_{shuffleIndex}_nCosines_{nCosines}"

    if not (len(qmTotalEnergy) == len(qmTorsionEnergy) == len(mmTotalEnergy) == len(mmFittedTorsionEnergy)):
        raise ValueError("QM/MM arrays must have matching lengths.")

    angleStepDegrees = 10
    angles = np.arange(len(qmTotalEnergy), dtype=float) * angleStepDegrees

    # Increased width ratio for the table column (0.6 instead of 0.4) to accommodate bigger tables
    fig = plt.figure(figsize=(18, 9))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 0.75])
    
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    ax_table = fig.add_subplot(gs[2])
    ax_table.axis('off')

    # --- Plotting Logic (Same as before) ---
    for n in range(1, 7):
        ax0.plot(angles, qmTotalEnergy, color=magenta, linewidth=1+n, alpha=0.1)
        ax0.plot(angles, mmTotalEnergy, color=brightCyan, linewidth=1+n, alpha=0.1)
    ax0.plot(angles, qmTotalEnergy, label='QM target', linewidth=2, color=magenta)
    ax0.plot(angles, mmTotalEnergy, label='MM fit', linewidth=2, color=brightCyan)
    ax0.set_title("Total Energies", fontsize=14)
    ax0.set_xlabel('Torsion Angle')
    ax0.set_ylabel('Energy (Kcal / mol)')
    ax0.set_xlim(0, 360)
    ax0.set_xticks(np.arange(0, 361, 60))
    ax0.legend(loc="upper right")

    for n in range(1, 7):
        ax1.plot(angles, qmTorsionEnergy, color=magenta, linewidth=1+n, alpha=0.1)
        ax1.plot(angles, mmFittedTorsionEnergy, color=brightCyan, linewidth=1+n, alpha=0.1)
    ax1.plot(angles, qmTorsionEnergy, label='QM target', linewidth=2, color=magenta)
    ax1.plot(angles, mmFittedTorsionEnergy, label='MM fit', linewidth=2, color=brightCyan)
    
    if isinstance(cosineComponents, dict):
        componentItems = cosineComponents.items()
    else:
        componentItems = cosineComponents
    for freq, cosineComponent in componentItems:
        ax1.plot(angles, cosineComponent, linestyle='--', color=white, alpha=0.4)
    
    ax1.plot([], [], linestyle='--', color=white, alpha=0.5, label='Parameters')
    ax1.set_title("Torsion Energies", fontsize=14)
    ax1.set_xlabel('Torsion Angle')
    ax1.set_xlim(0, 360)
    ax1.set_xticks(np.arange(0, 361, 60))
    ax1.legend(loc="upper right")

    # --- Table Logic ---
    
    # Display Tolerance at the very top
    ax_table.text(0.5, 0.98, f"Target Tolerance: {tol:.4f}", 
                  color="cyan", fontsize=14, fontweight='bold', ha='center', va='top', family='monospace')

    def create_table_data(metrics, mae):
        return [
            ["Score", f"{metrics['composite_score']:.3f}"],
            ["Loc", f"{metrics['location_score']:.3f}"],
            ["Amp", f"{metrics['amplitude_score']:.3f}"],
            ["Count", f"{metrics['stationary_count_score']:.3f}"],
            ["nMAE", f"{metrics['normalized_mae_score']:.3f}"],
            ["MAE", f"{mae:.3f}"]
        ]

    # Helper to style the table according to your specific rules
    def style_metric_table(table, tolerance):
        table.auto_set_font_size(False)
        table.set_fontsize(14) # Bigger text
        
        for (row, col), cell in table.get_celld().items():
            # Background and Edges
            cell.set_facecolor("#151515")
            cell.set_edgecolor("#444444")
            
            # Text properties
            cell_text = cell.get_text()
            cell_text.set_family("monospace")
            
            # Rule 1: Text (Labels in Col 0) are Yellow
            if col == 0:
                cell_text.set_color("yellow")
                # Rule 2: Score heading (Row 0) is Bold
                if row == 0:
                    cell_text.set_weight("bold")
            elif row == 5:
                cell_text.set_color("white")
                cell_text.set_weight("bold")
            
            # Rule 3: Numeric Values (Col 1)
            else:
                try:
                    val = float(cell_text.get_text())
                    # Rule 4: Green < tol, Orange < 2*tol, Red > 2*tol
                    if val < tolerance:
                        color = "lime"
                    elif val < (2 * tolerance):
                        color = "orange"
                    else:
                        color = "red"
                    
                    cell_text.set_color(color)
                    # Rule 2: Score value (Row 0, Col 1) also Bold
                    if row == 0:
                        cell_text.set_weight("bold")
                except ValueError:
                    cell_text.set_color("white")

    # Position Table 1: Total Metrics
    ax_table.text(0.1, 0.90, "TOTAL METRICS", color=white, fontweight='bold', ha='left', fontsize=12)
    table_total = ax_table.table(
        cellText=create_table_data(totalMetrics, maeTotal),
        colWidths=[0.5, 0.4],
        cellLoc='center',
        bbox=[0.1, 0.60, 0.8, 0.28] # Taller bbox for "bigger" tables
    )
    style_metric_table(table_total, tol)

    # Position Table 2: Torsion Metrics
    ax_table.text(0.1, 0.52, "TORSION METRICS", color=white, fontweight='bold', ha='left', fontsize=12)
    table_torsion = ax_table.table(
        cellText=create_table_data(torsionMetrics, maeTorsion),
        colWidths=[0.5, 0.4],
        cellLoc='center',
        bbox=[0.1, 0.22, 0.8, 0.28]
    )
    style_metric_table(table_torsion, tol)

    # Convergence status
    conv_color = "lime" if converged else "red"
    conv_text = "CONVERGED" if converged else "NOT CONVERGED"
    ax_table.text(
        0.5, 0.08, conv_text,
        color=conv_color, fontsize=18, fontweight="bold", ha="center",
        bbox=dict(facecolor="none", edgecolor=conv_color, boxstyle="round,pad=0.5")
    )

    plt.tight_layout()
    # Save with figure facecolor to maintain dark theme if set_rc_params defines it
    fig.savefig(p.join(outDir, f"{shuffleLabel}.png"), facecolor=fig.get_facecolor())
    plt.close(fig)


#####################################################################

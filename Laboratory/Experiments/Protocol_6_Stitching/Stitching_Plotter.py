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
from . import Stitching_Assistant

## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass

def plot_run_mean_average_error(outDir: DirectoryPath, fittingScoresCsv: FilePath) -> None:
    """Plot the run-wide mean average error curves."""

    maeDf = pd.read_csv(fittingScoresCsv)

    set_rc_params()
    ## init some colors to be used
    white: str = '#FFFFFF'
    brightGreen: str = '#00FF00'
    brightOrange: str = '#FFA500'
    fig, axis = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))  # create figure


    allTorsionsDf = maeDf[maeDf["torsion_tag"] == "All_Torsions"]
    individualTorsionDf = maeDf[maeDf["torsion_tag"] != "All_Torsions"]

    # Determine which columns to plot (backwards compatible)
    torsion_score_col = "fit_score_torsion" if "fit_score_torsion" in allTorsionsDf.columns else "mae_torsion"
    total_score_col = "fit_score_total" if "fit_score_total" in allTorsionsDf.columns else "mae_total"

    # ==================== COMPOSITE SCORES (opaque glowing) ====================
    x = allTorsionsDf["shuffle"]
    # glow background
    for n in range(1, 7):
        axis.plot(x, allTorsionsDf[torsion_score_col], color=brightGreen, linewidth=1 + n * 1.25, alpha=0.12)
    axis.plot(x, allTorsionsDf[torsion_score_col], label='Torsion Score', linewidth=2.5, color=brightGreen, alpha=1.0)

    for n in range(1, 7):
        axis.plot(x, allTorsionsDf[total_score_col], color=brightOrange, linewidth=1 + n * 1.25, alpha=0.12)
    axis.plot(x, allTorsionsDf[total_score_col], label='Total Score', linewidth=2.5, color=brightOrange, alpha=1.0)

    # subtle traces for individual torsions (composite score)
    for torsionTag, df in individualTorsionDf.groupby("torsion_tag"):
        col = "fit_score_torsion" if "fit_score_torsion" in df.columns else "mae_torsion"
        axis.plot(df["shuffle"], df[col], linewidth=1, color=brightGreen, alpha=0.35)
        col_total = "fit_score_total" if "fit_score_total" in df.columns else "mae_total"
        axis.plot(df["shuffle"], df[col_total], linewidth=1, color=brightOrange, alpha=0.25)

    # ==================== SUBCOMPONENTS (thin, semi-transparent) ====================
    comp_cols = [
        ("torsion_location_score", "Location"),
        ("torsion_amplitude_score", "Amplitude"),
        ("torsion_stationary_count_score", "Count"),
        ("torsion_normalized_mae_score", "nMAE"),
    ]
    comp_colors = ["#FFD700", "#FF69B4", "#87CEFA", "#ADFF2F"]
    for (col, label), ccol in zip(comp_cols, comp_colors):
        if col in allTorsionsDf.columns:
            axis.plot(x, allTorsionsDf[col], linewidth=1.0, color=ccol, alpha=0.6, linestyle='-')

    comp_cols_total = [
        ("total_location_score", "Total_Loc"),
        ("total_amplitude_score", "Total_Amp"),
        ("total_stationary_count_score", "Total_Count"),
        ("total_normalized_mae_score", "Total_nMAE"),
    ]
    for (col, label), ccol in zip(comp_cols_total, comp_colors):
        if col in allTorsionsDf.columns:
            axis.plot(x, allTorsionsDf[col], linewidth=1.0, color=ccol, alpha=0.4, linestyle='--')

    # Add labels and legend
    axis.set_xlabel('Number of Fitting Shuffles')
    axis.set_ylabel('Composite Fit Score (lower is better)')
    axis.set_title('Fit Scores vs Number of Fitting Shuffles')

    torsionFlatLined, totalFlatLined = Stitching_Assistant.score_flatline_status(maeDf, torsionTag="All_Torsions")
    if torsionFlatLined and totalFlatLined:
        axis.text(
            0.5, 0.92,
            "FLATLINED (TORSION + TOTAL)",
            transform=axis.transAxes,
            color="red",
            fontsize=14,
            fontweight="bold",
            ha="center",
            bbox=dict(facecolor="none", edgecolor="red", boxstyle="round,pad=0.4"),
        )

    # Build a concise legend: include main composite traces
    handles, labels = axis.get_legend_handles_labels()
    filtered_handles_labels = [(h, l) for h, l in zip(handles, labels) if l in ['Torsion Score', 'Total Score']]
    if filtered_handles_labels:
        handles, labels = zip(*filtered_handles_labels)
        axis.legend(handles, labels)
    else:
        axis.legend()

    # Get unique shuffle values and ensure they are integers
    if not allTorsionsDf.empty:
        shuffle_values = allTorsionsDf["shuffle"].unique().astype(int)
        axis.set_xticks(shuffle_values)  # Set all integer ticks initially
        axis.xaxis.set_major_locator(MaxNLocator(integer=True, nbins='auto'))  # Auto-adjust ticks
        current_ticks = axis.get_xticks()
        valid_ticks = current_ticks[(current_ticks >= 0) & (current_ticks <= shuffle_values[-1])]
        if shuffle_values[-1] not in valid_ticks:
            valid_ticks = np.append(valid_ticks, shuffle_values[-1])
        axis.set_xticks(np.sort(valid_ticks))

    plt.savefig(p.join(outDir, "run_mean_average_error.png"), bbox_inches='tight')
    plt.close()



def plot_mean_average_error(torsionFittingDir: DirectoryPath, fittingScoresCsv: FilePath, torsionTag: str) -> None:
    """Plot the mean average error curves for one torsion."""

    ## get data from csv
    maeDf = pd.read_csv(fittingScoresCsv)
    torsionMaeDf = maeDf[maeDf["torsion_tag"] == torsionTag]
    

    set_rc_params()
    ## init some colors to be used
    brightGreen: str = '#00FF00'
    brightOrange: str = '#FFA500'

    fig, axis = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))  # create figure

    # Choose columns with fallbacks
    torsion_score_col = "fit_score_torsion" if "fit_score_torsion" in torsionMaeDf.columns else "mae_torsion"
    total_score_col = "fit_score_total" if "fit_score_total" in torsionMaeDf.columns else "mae_total"

    x = torsionMaeDf["shuffle"] if "shuffle" in torsionMaeDf.columns else np.arange(len(torsionMaeDf))

    ## plot torsion composite score with GLOW (opaque)
    for n in range(1, 7):
        axis.plot(x, torsionMaeDf[torsion_score_col], color=brightGreen, linewidth=1 + n * 1.25, alpha=0.12)
    axis.plot(x, torsionMaeDf[torsion_score_col], label='Torsion Score', linewidth=2.5, color=brightGreen, alpha=1.0)

    ## plot total composite score with GLOW (opaque)
    for n in range(1, 7):
        axis.plot(x, torsionMaeDf[total_score_col], color=brightOrange, linewidth=1 + n * 1.25, alpha=0.12)
    axis.plot(x, torsionMaeDf[total_score_col], label='Total Score', linewidth=2.5, color=brightOrange, alpha=1.0)

    # Plot subcomponents for this torsion (thin semi-transparent lines)
    comp_cols = [
        ("torsion_location_score", "Location"),
        ("torsion_amplitude_score", "Amplitude"),
        ("torsion_stationary_count_score", "Count"),
        ("torsion_normalized_mae_score", "nMAE"),
    ]
    comp_colors = ["#FFD700", "#FF69B4", "#87CEFA", "#ADFF2F"]
    for (col, label), ccol in zip(comp_cols, comp_colors):
        if col in torsionMaeDf.columns:
            axis.plot(x, torsionMaeDf[col], linewidth=1.0, color=ccol, alpha=0.6)

    # Also show the total subcomponents (dashed, fainter)
    comp_cols_total = [
        ("total_location_score", "Total_Loc"),
        ("total_amplitude_score", "Total_Amp"),
        ("total_stationary_count_score", "Total_Count"),
        ("total_normalized_mae_score", "Total_nMAE"),
    ]
    for (col, label), ccol in zip(comp_cols_total, comp_colors):
        if col in torsionMaeDf.columns:
            axis.plot(x, torsionMaeDf[col], linewidth=1.0, color=ccol, alpha=0.4, linestyle='--')

    # Add labels and legend
    axis.set_xlabel('Number of Fitting Shuffles')
    axis.set_ylabel('Composite Fit Score (lower is better)')
    axis.set_title('Fit Scores vs Number of Fitting Shuffles')

    torsionFlatLined, totalFlatLined = Stitching_Assistant.score_flatline_status(maeDf, torsionTag=torsionTag)
    if torsionFlatLined and totalFlatLined:
        axis.text(
            0.5, 0.92,
            "FLATLINED (TORSION + TOTAL)",
            transform=axis.transAxes,
            color="red",
            fontsize=14,
            fontweight="bold",
            ha="center",
            bbox=dict(facecolor="none", edgecolor="red", boxstyle="round,pad=0.4"),
        )
    axis.legend()

    plt.savefig(p.join(torsionFittingDir, "mean_average_error.png"), bbox_inches='tight')
    plt.close()
#####################################################################
def _extract_number(filename: str) -> int:
    """Extract the numeric suffix from a filename."""
    # Extract the number from the filename
    return int(filename.split('_')[-1].split('.')[0])


def _extract_fit_frame_key(filename: str) -> tuple[int, int]:
    """Sort fitting frames by cosine count and shuffle index."""
    match = re.search(r"fitting_shuffle_(\d+)_nCosines_(\d+)\.png$", os.path.basename(filename))
    if match:
        return int(match.group(2)), int(match.group(1))
    return (0, _extract_number(filename))
    
#####################################################################

def make_gif(inDir: DirectoryPath, outGif: FilePath, batchSize: int = 50, duration: int = 100) -> None:
    """Create a GIF from fitting PNG frames."""
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
def set_rc_params() -> None:
    """Set matplotlib defaults for stitching plots."""
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
def plot_qmmm_energies(qmTotalEnergy: np.ndarray,
                        qmTorsionEnergy: np.ndarray,
                        mmTotalEnergy: np.ndarray,
                        mmFittedTorsionEnergy: np.ndarray,
                        cosineComponents: dict,
                        torsionMetrics: dict,
                        totalMetrics: dict,
                        maeTorsion: float,
                        maeTotal: float,
                        nCosines: int,
                        converged: bool,
                        outDir: DirectoryPath,
                        shuffleIndex: int,
                        tol: float,
                        qmTotalEnergyRaw: np.ndarray | None = None,
                        mmTotalEnergyRaw: np.ndarray | None = None) -> None:
    """Plot QM, MM, and fitted torsion energies for one shuffle."""
    
    set_rc_params()
    white :str = '#FFFFFF'
    magenta : str = '#FF00FF'
    brightCyan: str = '#00FFFF'
    shuffleLabel = f"fitting_shuffle_{shuffleIndex}_nCosines_{nCosines}"

    if not (len(qmTotalEnergy) == len(qmTorsionEnergy) == len(mmTotalEnergy) == len(mmFittedTorsionEnergy)):
        raise ValueError("QM/MM arrays must have matching lengths.")
    if qmTotalEnergyRaw is not None and len(qmTotalEnergyRaw) != len(qmTotalEnergy):
        raise ValueError("qmTotalEnergyRaw must match qmTotalEnergy length.")
    if mmTotalEnergyRaw is not None and len(mmTotalEnergyRaw) != len(mmTotalEnergy):
        raise ValueError("mmTotalEnergyRaw must match mmTotalEnergy length.")

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
    if qmTotalEnergyRaw is not None:
        ax0.plot(
            angles,
            qmTotalEnergyRaw,
            linestyle=':',
            linewidth=1.0,
            color=magenta,
            alpha=0.3,
            label='QM raw (unsmoothed)',
        )
    if mmTotalEnergyRaw is not None:
        ax0.plot(
            angles,
            mmTotalEnergyRaw,
            linestyle=':',
            linewidth=1.0,
            color=brightCyan,
            alpha=0.3,
            label='MM raw (unsmoothed)',
        )
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

    def create_table_data(metrics: dict, mae: float) -> list[list[str]]:
        return [
            ["Score", f"{metrics['composite_score']:.3f}"],
            ["Loc", f"{metrics['location_score']:.3f}"],
            ["Amp", f"{metrics['amplitude_score']:.3f}"],
            ["Count", f"{metrics['stationary_count_score']:.3f}"],
            ["nMAE", f"{metrics['normalized_mae_score']:.3f}"],
            ["MAE", f"{mae:.3f}"]
        ]

    # Helper to style the table according to your specific rules
    def style_metric_table(table, tolerance: float) -> object:
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

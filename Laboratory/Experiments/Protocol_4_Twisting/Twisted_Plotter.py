import os
from os import path as p
import pandas as pd
import numpy as np
## PLOTTING LIBRARIES ##
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from typing import List, Optional # Added for type hinting

# Placeholder classes (extend if needed)
class FilePath:
    pass
class DirectoryPath:
    pass

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def twist_plotting_protocol(scan_df: List[pd.DataFrame], 
                             scan_averages_df: pd.DataFrame, 
                             sp_df: Optional[List[pd.DataFrame]], 
                             sp_averages_df: Optional[pd.DataFrame], 
                             out_dir: DirectoryPath, 
                             torsion_tag: str, 
                             config: dict) -> dict:
    """
    Main protocol for generating plots related to torsion scanning.

    Args:
        scan_df: List of DataFrames, each containing individual scan energies.
        scan_averages_df: DataFrame with averaged scan energies.
        sp_df: Optional list of DataFrames, each with individual single point energies.
        sp_averages_df: Optional DataFrame with averaged single point energies.
        out_dir: The output directory where plots will be saved.
        torsion_tag: Identifier for the torsion being plotted.
        config: Configuration dictionary.

    Returns:
        The updated configuration dictionary with paths to generated plots.
    """
    ## make a dir to store plots
    scanPlotDir = p.join(out_dir, "plots") # type: ignore
    os.makedirs(scanPlotDir, exist_ok=True)

    set_rc_params()

    ## extract method for scan 
    scanMethod = config["torsionScanInfo"]["scanMethod"]
    scanSolvation = config["torsionScanInfo"]["scanSolvationMethod"]
    if scanSolvation is None:
        scanQmMethod = scanMethod
    else:
        scanQmMethod = scanMethod + " [ " + scanSolvation + " ]"

    scanPng = plot_individual_vs_average(scan_df, scan_averages_df, scanPlotDir, torsion_tag, "scan", scanQmMethod)
    config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"][torsion_tag]["scanPng"] = scanPng


    if sp_df is None or sp_averages_df is None or config["torsionScanInfo"]["singlePointMethod"] is None : # Added sp_df/sp_averages_df check
        config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"][torsion_tag]["spPng"] = None
        config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"][torsion_tag]["scanVsSpPng"] = None

    else:
        ## extract qm method for single point
        singlepointMethod = config["torsionScanInfo"]["singlePointMethod"]
        singlepointSolvation = config["torsionScanInfo"]["singlePointSolvationMethod"]
        if singlepointSolvation is None:
            singlepointQmMethod = singlepointMethod
        else:
            singlepointQmMethod = singlepointMethod + " [ " + singlepointSolvation + " ]"

        spPng = plot_individual_vs_average(sp_df, sp_averages_df, scanPlotDir, torsion_tag, "SP", singlepointQmMethod)
        scanVsSpPng = plot_scan_singlepoint_comparison(scan_averages_df, sp_averages_df, scanPlotDir, torsion_tag, scanQmMethod, singlepointQmMethod)
        config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"][torsion_tag]["spPng"] = spPng
        config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"][torsion_tag]["scanVsSpPng"] = scanVsSpPng

    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def set_rc_params() -> None:
    """Sets matplotlib rcParams for consistent plot styling."""
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


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_individual_vs_average(individual_dfs: List[pd.DataFrame], 
                                average_df: pd.DataFrame, 
                                out_dir: DirectoryPath, 
                                torsion_tag: str, 
                                tag: str, 
                                qm_method: str) -> FilePath:
    """
    Plots individual scan/SP energies against their average.

    Args:
        individual_dfs: List of DataFrames, each with 'Angle' and 'Energy' columns.
        average_df: DataFrame with 'Angle' and a column named by torsion_tag for average energies.
        out_dir: Directory to save the plot.
        torsion_tag: Identifier for the torsion.
        tag: Type of plot ('scan' or 'SP').
        qm_method: QM method string for the plot title.

    Returns:
        Path to the generated plot PNG file.
    """
    ## init some colors to be used
    white :str = '#FFFFFF'
    brightGreen: str = '#00FF00'
    magenta : str = '#FF00FF'

    if tag == "scan":
        averageLineColor = brightGreen
        saveTag = "01_scan_energies"
    elif tag == "SP":
        averageLineColor = magenta
        saveTag = "02_SP_energies"
    else: # Default or error case
        averageLineColor = white
        saveTag = "unknown_energies"


    plt.figure(figsize=(12, 8))
    for df in individual_dfs:
        plt.plot(df['Angle'], df["Energy"], color=white, alpha=0.3)

    for n in range(1, 7):
        plt.plot(average_df['Angle'], average_df[torsion_tag], label='Average', color=averageLineColor, linewidth=1+n, alpha=0.1)
    plt.plot(average_df['Angle'], average_df[torsion_tag], label='Average', color=averageLineColor, linewidth=1)

    plt.xlabel('Angle (degrees)')
    plt.ylabel('Energy (Kcal / mol)')
    plt.title(f'Torsion scans for {torsion_tag} at {qm_method}')


    plotPng = p.join(out_dir,f"{saveTag}_{torsion_tag}.png") # type: ignore
    # Save the plot to the specified directory
    plt.savefig(plotPng)
    plt.close()

    return plotPng # type: ignore


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_scan_singlepoint_comparison(scan_df: pd.DataFrame, 
                                     sp_df: pd.DataFrame, 
                                     out_dir: DirectoryPath, 
                                     torsion_tag: str, 
                                     scan_qm_method: str, 
                                     singlepoint_qm_method: str) -> FilePath:
    """
    Plots a comparison of averaged scan energies vs averaged single point energies.

    Args:
        scan_df: DataFrame with averaged scan energies ('Angle', torsion_tag column).
        sp_df: DataFrame with averaged single point energies ('Angle', torsion_tag column).
        out_dir: Directory to save the plot.
        torsion_tag: Identifier for the torsion.
        scan_qm_method: QM method string for scan data.
        singlepoint_qm_method: QM method string for single point data.

    Returns:
        Path to the generated plot PNG file.
    """
    ## init some colors to be used
    white :str = '#FFFFFF' # Unused, but kept from original
    brightGreen: str = '#00FF00'
    magenta : str = '#FF00FF'

    plt.figure(figsize=(12, 8))
    

    ## plot scan average with GLOW
    for n in range(1, 7):
        plt.plot(scan_df['Angle'], scan_df[torsion_tag], color=brightGreen, linewidth=1+n, alpha=0.1)
    plt.plot(scan_df['Angle'], scan_df[torsion_tag], 
             label='Scan Data', color=brightGreen)
    ## plot single point average with GLOW
    for n in range(1, 7):
        plt.plot(sp_df['Angle'], sp_df[torsion_tag], color=magenta, linewidth=1+n, alpha=0.1)
    plt.plot(sp_df['Angle'], sp_df[torsion_tag], 
             label='Single Point Data', color=magenta)
    
    # Add labels and legend
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Energy (Kcal / mol)')
    plt.title('Scan vs Single Point Energy Comparison')
    plt.legend(
        [f"Scan at {scan_qm_method}", f"Single Point at {singlepoint_qm_method}"],
        loc='best',   
        handlelength=0,
        handletextpad=0,
        labelcolor=[brightGreen, magenta]
    )
    
    plotPng = p.join(out_dir,f"03_scan_vs_SP_{torsion_tag}.png") # type: ignore
    # Save the plot to the specified directory
    plt.savefig(plotPng)
    plt.close()

    return plotPng

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

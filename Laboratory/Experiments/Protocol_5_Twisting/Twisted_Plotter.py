import os
from os import path as p
import pandas as pd
import numpy as np
## PLOTTING LIBRARIES ##
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def twist_plotting_protocol(scanDf, scanAveragesDf, spDf, spAveragesDf, outDir, torsionTag, config):
    ## make a dir to store plots
    scanPlotDir = p.join(outDir, "plots")
    os.makedirs(scanPlotDir, exist_ok=True)

    set_rc_params()

    ## extract method for scan 
    scanMethod = config["torsionScanInfo"]["scanMethod"]
    scanSolvation = config["torsionScanInfo"]["scanSolvationMethod"]
    if scanSolvation is None:
        scanQmMethod = scanMethod
    else:
        scanQmMethod = scanMethod + " [ " + scanSolvation + " ]"

    scanPng = plot_individual_vs_average(scanDf, scanAveragesDf, scanPlotDir, torsionTag, "scan", scanQmMethod)
    config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"][torsionTag]["scanPng"] = scanPng


    if  config["torsionScanInfo"]["singlePointMethod"] is None:
        config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"][torsionTag]["spPng"] = None
        config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"][torsionTag]["scanVsSpPng"] = None

    else:
        ## extract qm method for single point
        singlepointMethod = config["torsionScanInfo"]["singlePointMethod"]
        singlepointSolvation = config["torsionScanInfo"]["singlePointSolvationMethod"]
        if singlepointSolvation is None:
            singlepointQmMethod = singlepointMethod
        else:
            singlepointQmMethod = singlepointMethod + " [ " + singlepointSolvation + " ]"

        spPng = plot_individual_vs_average(spDf, spAveragesDf, scanPlotDir, torsionTag, "SP", singlepointQmMethod)
        scanVsSpPng = plot_scan_singlepoint_comparison(scanAveragesDf, spAveragesDf, scanPlotDir, torsionTag, scanQmMethod, singlepointQmMethod)
        config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"][torsionTag]["spPng"] = spPng
        config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"][torsionTag]["scanVsSpPng"] = scanVsSpPng

    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

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


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_individual_vs_average(individualDfs, averageDf, outDir, torsionTag, tag, qmMethod):

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

    plt.figure(figsize=(12, 8))
    for df in individualDfs:
        plt.plot(df['Angle'], df["Energy"], color=white, alpha=0.3)

    for n in range(1, 7):
        plt.plot(averageDf['Angle'], averageDf[torsionTag], label='Average', color=averageLineColor, linewidth=1+n, alpha=0.1)
    plt.plot(averageDf['Angle'], averageDf[torsionTag], label='Average', color=averageLineColor, linewidth=1)

    plt.xlabel('Angle (degrees)')
    plt.ylabel('Energy (Kcal / mol)')
    plt.title(f'Torsion scans for {torsionTag} at {qmMethod}')


    plotPng = p.join(outDir,f"{saveTag}_{torsionTag}.png")
    # Save the plot to the specified directory
    plt.savefig(plotPng)
    plt.close()

    return plotPng


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_scan_singlepoint_comparison(scanDf, spDf, outDir, torsionTag, scanQmMethod, singlepointQmMethod):

    ## init some colors to be used
    white :str = '#FFFFFF'
    brightGreen: str = '#00FF00'
    magenta : str = '#FF00FF'

    plt.figure(figsize=(12, 8))
    

    ## plot scan average with GLOW
    for n in range(1, 7):
        plt.plot(scanDf['Angle'], scanDf[torsionTag], color=brightGreen, linewidth=1+n, alpha=0.1)
    plt.plot(scanDf['Angle'], scanDf[torsionTag], 
             label='Scan Data', color=brightGreen)
    ## plot single point average with GLOW
    for n in range(1, 7):
        plt.plot(spDf['Angle'], spDf[torsionTag], color=magenta, linewidth=1+n, alpha=0.1)
    plt.plot(spDf['Angle'], spDf[torsionTag], 
             label='Single Point Data', color=magenta)
    
    # Add labels and legend
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Energy (Kcal / mol)')
    plt.title('Scan vs Single Point Energy Comparison')
    plt.legend(
        [f"Scan at {scanQmMethod}", f"Single Point at {singlepointQmMethod}"],
        loc='best',   
        handlelength=0,
        handletextpad=0,
        labelcolor=[brightGreen, magenta]
    )
    
    plotPng = p.join(outDir,f"03_scan_vs_SP_{torsionTag}.png")
    # Save the plot to the specified directory
    plt.savefig(plotPng)
    plt.close()

    return plotPng

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

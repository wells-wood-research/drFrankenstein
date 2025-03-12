import os
from os import path as p
import pandas as pd
import numpy as np
## PLOTTING LIBRARIES ##
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def _debug_plot_array(data, outDir):
    plt.figure(figsize=(12, 8))
    plt.plot(data)
    plt.savefig(p.join(outDir, "debug.png"))
    plt.close()

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_scan_singlepoint_comparison(scanDf, spDf, outDir, tag):

    plt.figure(figsize=(12, 8))
    
    # Plot the first DataFrame
    plt.plot(scanDf['Angle'], scanDf[tag], 
             label='Scan Data', color='blue')
    
    # Plot the second DataFrame
    plt.plot(spDf['Angle'], spDf[tag], 
             label='Single Point Data', color='red')
    
    # Add labels and legend
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Energy')
    plt.title('Scan vs Single Point Energy Comparison')
    plt.legend()
    
    # Save the plot to the specified directory
    plt.savefig(p.join(outDir,f"singlepoint_comparison_{tag}.png"))
    plt.close()
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲


def plot_torsion_scans(torsionTopDir, scanRawDfs, scanAverageDf, rollingAverageDf, meanAverageErrors):
    plotDir = p.join(torsionTopDir, "plots")
    os.makedirs(plotDir, exist_ok=True)

    plot_raw_data(scanRawDfs, plotDir)
    plot_average_data(scanAverageDf, plotDir)
    plot_rolling_average_data(rollingAverageDf, plotDir)
    plot_mean_average_data(meanAverageErrors, plotDir)

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_raw_data(scanRawDfs, plotDir):
    colors = cm.YlGn(np.linspace(0, 1, len(scanRawDfs)))
    
    plt.figure(figsize=(12, 8))
    
    for (batchIndex, df), color in zip(scanRawDfs.items(), colors):
        for column in df.columns:
            if column != 'Angle':
                plt.scatter(df['Angle'], df[column], label=f'{batchIndex} - {column}', color=color)
    
    plt.xticks(np.arange(-180, 181, 60))
    plt.xlabel('Angle')
    plt.ylabel('Relative Energy (kcal/mol)')
    plt.title('Raw Data')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(p.join(plotDir, "raw_data.png"))
    plt.close()

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_average_data(scanAverageDf, plotDir):
    colors = cm.YlGn(np.linspace(0, 1, len(scanAverageDf.columns)-1))

    plt.figure(figsize=(12, 8))
    
    for column, color in zip(scanAverageDf.columns, colors):
        if column != 'Angle':
            plt.plot(scanAverageDf['Angle'], scanAverageDf[column], label=column, color=color)

    plt.xticks(np.arange(-180, 181, 60))
    plt.xlabel('Angle')
    plt.ylabel('Relative Energy (kcal/mol)')
    plt.title('Average Data')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(p.join(plotDir, "average_data.png"))
    plt.close()
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_rolling_average_data(rollingAverageDf, plotDir):
    colors = cm.YlGn(np.linspace(0, 1, len(rollingAverageDf.columns)-1))

    plt.figure(figsize=(12, 8))
    
    for column, color in zip(rollingAverageDf.columns, colors):
        if column != 'Angle':
            plt.plot(rollingAverageDf['Angle'], rollingAverageDf[column], label=column, color=color)

    plt.xticks(np.arange(-180, 181, 60))

    plt.xlabel('Angle')
    plt.ylabel('Relative Energy (kcal/mol)')
    plt.title('Rolling Average Data')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(p.join(plotDir, "rolling_average_data.png"))
    plt.close()
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_mean_average_data(meanAverageErrors, plotDir):
    plt.figure(figsize=(12, 8))
    
    plt.plot(meanAverageErrors)
    plt.ylabel('Relative Energy (kcal/mol)')
    plt.title('Mean Average Error')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(p.join(plotDir, "mean_average_error.png"))
    plt.close()

    #🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

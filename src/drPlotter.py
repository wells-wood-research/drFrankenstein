import os
from os import path as p
import pandas as pd
import numpy as np
## PLOTTING LIBRARIES ##
import matplotlib.pyplot as plt
import matplotlib.cm as cm


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

import os
from os import path as p
import pandas as pd
from functools import reduce
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text

def inputs():
    benchmarkDir = '/home/esp/scriptDevelopment/drFrankenstein/Benchmarking'
    return benchmarkDir

def main():
    benchmarkDir = inputs()
    chargesDf = get_data(benchmarkDir)
    outPng = p.join(benchmarkDir, 'violins.png')
    plot_swarm(chargesDf, outPng)

def get_data(benchmarkDir):
    allDataDf = pd.DataFrame(columns=[])
    for methodName in os.listdir(benchmarkDir):
        methodDir = p.join(benchmarkDir, methodName)
        if not p.isdir(methodDir):
            continue
        chargeFittingDir = p.join(methodDir, 'charge_calculations', 'charge_fitting')
        chargesCsv = p.join(chargeFittingDir, 'charges.csv')
        chargesDf = pd.read_csv(chargesCsv, index_col='Unnamed: 0')
        allDataDf[methodName] = chargesDf['Charge']
    allDataDf['atomIndex'] = chargesDf['atomIndex']
    allDataDf['atomElement'] = chargesDf['atomElement']
    return allDataDf

def plot_swarm(df, outPng):
    meltedDf = df.melt(id_vars=['atomIndex', 'atomElement'], value_vars=['XTB2-revPBE-SVP', 'XTB2-HF', 'XTB2-revPBE-TZVP', 'HF-HF'], var_name='Method', value_name='Value')
    markers = {'C': 's', 'N': '^', 'O': 'o', 'H': 'o'}
    colors = {'C': 'black', 'N': 'blue', 'O': 'red', 'H': 'green'}
    plt.figure(figsize=(12, 6))
    for element in colors.keys():
        subset = meltedDf[meltedDf['atomElement'] == element]
        sns.swarmplot(x='Method', y='Value', data=subset, color=colors[element], marker=markers[element], label=element, dodge=True, size=10)
    (handles, labels) = plt.gca().get_legend_handles_labels()
    byLabel = dict(zip(labels, handles))
    plt.legend(byLabel.values(), byLabel.keys(), title='Atom Element')
    plt.title('Swarm Plot of Methods')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(outPng)
    plt.close()
main()
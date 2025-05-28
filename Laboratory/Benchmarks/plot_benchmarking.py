import os
from os import path as p
import pandas as pd
from functools import reduce
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text


def inputs() -> str:
    """Returns the benchmark directory path."""
    benchmarkDir: str = "/home/esp/scriptDevelopment/drFrankenstein/Benchmarking"

    return benchmarkDir

def main() -> None:
    """Main function to generate and save the swarm plot."""
    benchmarkDir: str = inputs()
    chargesDf: pd.DataFrame = get_data(benchmarkDir)
    outPng: str = p.join(benchmarkDir, "violins.png")
    plot_swarm(chargesDf, outPng)

def get_data(benchmarkDir: str) -> pd.DataFrame:
    """
    Reads charge data from CSV files in different method directories and compiles them into a single DataFrame.

    Args:
        benchmarkDir: The root directory containing subdirectories for each method.

    Returns:
        A pandas DataFrame containing the combined charge data.
    """
    allDataDf: pd.DataFrame = pd.DataFrame(columns=[])
    for methodName in os.listdir(benchmarkDir):
        methodDir: str = p.join(benchmarkDir, methodName)
        if not p.isdir(methodDir):
            continue
        chargeFittingDir: str = p.join(methodDir, "charge_calculations", "charge_fitting")
        chargesCsv: str = p.join(chargeFittingDir, "charges.csv")
        chargesDf: pd.DataFrame = pd.read_csv(chargesCsv, index_col="Unnamed: 0")
        allDataDf[methodName] = chargesDf["Charge"]

    allDataDf["atomIndex"]  = chargesDf["atomIndex"]
    allDataDf["atomElement"]  = chargesDf["atomElement"]

    return allDataDf

def plot_swarm(df: pd.DataFrame, outPng: str) -> None:
    """
    Generates a swarm plot of charge values for different methods and atom types and saves it to a file.

    Args:
        df: DataFrame containing the charge data, with columns for atomIndex, atomElement, and each method.
        outPng: The path to save the output PNG image.
    """
    # Melt the DataFrame to long format for seaborn
    meltedDf: pd.DataFrame = df.melt(id_vars=['atomIndex', 'atomElement'], 
                        value_vars=['XTB2-revPBE-SVP', 'XTB2-HF', 
                                    'XTB2-revPBE-TZVP', 'HF-HF'],
                        var_name='Method', value_name='Value')
    markers: dict[str, str] = {'C': 's', 'N': '^', 'O': 'o', 'H': 'o'}
    colors: dict[str, str] = {'C': 'black', 'N': 'blue', 'O': 'red', 'H': 'green'}
    # Create the swarm plot with different markers
    plt.figure(figsize=(12, 6))
    # Plot each atom element with its specific marker and color
    for element in colors.keys():
        subset: pd.DataFrame = meltedDf[meltedDf['atomElement'] == element]
        sns.swarmplot(x='Method', y='Value', data=subset,
                       color=colors[element], marker=markers[element],
                      label=element, dodge=True, size = 10)
    # Customize the legend to show only unique values
    handles, labels = plt.gca().get_legend_handles_labels()
    byLabel: dict = dict(zip(labels, handles))
    plt.legend(byLabel.values(), byLabel.keys(), title='Atom Element')


    # Customize the plot
    plt.title('Swarm Plot of Methods')
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(outPng)
    plt.close()

main()
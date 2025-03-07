import numpy as np

import matplotlib
matplotlib.use('Agg')  # Set the Agg backend before importing pyplot
import matplotlib.pyplot as plt
import pandas as pd
import os
from os import path as p
from typing import Tuple, List
import plotly.graph_objects as go

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
## dummy classes
class FilePath:
    pass
class DirPath:
    pass

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def fourier_transform_protocol(qmTorsionEnergy, torsionTag, torsionFittingDir, sampleSpacing=10, maxFunctions=10):
    energyDataPadded: np.array = pad_energy_data(qmTorsionEnergy, paddingFactor=3)
    ## calculate signal length
    signalLength: int = len(energyDataPadded)
    ## run reverse fourier transform
    fftResult: np.array = perform_rfft(energyDataPadded)
    ## get frequencies, amplitudes and phases
    frequencies: np.array = get_frequencies(signalLength, sampleSpacing)
    amplitudes, phases = compute_amplitude_and_phase(fftResult, signalLength)
    ## construct angle x-axis
    angle: np.array = np.arange(signalLength) * sampleSpacing
    ## convert data to dataframe
    fourierDf: pd.DataFrame = convert_fourier_params_to_df(frequencies, amplitudes, phases)

    amberParamDf: pd.DataFrame = convert_params_to_amber_format(fourierDf)
    ## construct cosine components from parameters
    reconstructedSignal, cosineComponents, nFunctionsUsed = construct_cosine_components(amberParamDf, angle, maxFunctions)
    ## plot signal and components
    plot_signal_and_components(angle, qmTorsionEnergy, reconstructedSignal, cosineComponents, torsionFittingDir , torsionTag)
    ## write data to csv file
    outCsv: FilePath = p.join(torsionFittingDir, f"{torsionTag}.csv")
    amberParamDf.iloc[:nFunctionsUsed].to_csv(outCsv)

    return amberParamDf.iloc[:nFunctionsUsed]   
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def convert_fourier_params_to_df(frequencies: np.array, amplitudes: np.array, phases: np.array) -> pd.DataFrame:
    data = {"Frequency": frequencies, "Amplitude": amplitudes, "Phase": phases}
    dataDf = pd.DataFrame(data)
    dataDf = dataDf.sort_values(by='Amplitude', ascending=False).reset_index(drop=True)
    return dataDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def convert_params_to_amber_format(fourierDf: pd.DataFrame) -> pd.DataFrame:
    amberDf = pd.DataFrame()
    amberDf["Amplitude"] = fourierDf["Amplitude"]
    amberDf["Period"] = np.degrees(2 * np.pi * fourierDf["Frequency"])
    amberDf["Phase"] = np.degrees(fourierDf["Phase"]) * -1

    ## remove period == 0 (DC) Signal
    amberDf = amberDf[amberDf["Period"] > 0]

    amberDf.sort_values(by="Amplitude", ascending=False, inplace=True,ignore_index=True)
    return amberDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def construct_cosine_components(amberParamDf: pd.DataFrame,
                                 angle: np.array,
                                     maxFunctions: int,
                                       tolerance: float = 0.5) -> Tuple[np.array, List[Tuple[float, np.array]], int]:
    
    amplitudes = amberParamDf["Amplitude"]
    periods = amberParamDf["Period"]
    phases = amberParamDf["Phase"]

    sample_spacing = 10
    signalLength = 36
    
    # Angle array
    angle = np.arange(signalLength) * sample_spacing  # Shape: (N,)
    
    # Initialize reconstructed signal
    reconstructedSignal = np.zeros(signalLength)
    previousSignal = np.zeros(signalLength)

    ## init mean average error
    meanAverageError = np.inf
    ## collect cosine components for plotting
    cosineComponents = []
    # Construct each cosine component
    while True:
        for nFunctionsUsed in range(1, maxFunctions+1):
            amberComponent =  amplitudes[nFunctionsUsed] * (1 + np.cos(np.radians(periods[nFunctionsUsed] * angle - phases[nFunctionsUsed])))

            previousSignal = reconstructedSignal.copy()
            reconstructedSignal += amberComponent
            meanAverageError =  np.mean(np.abs(reconstructedSignal - previousSignal))
            cosineComponents.append((nFunctionsUsed , amberComponent))
            if meanAverageError < tolerance:
                break
        if meanAverageError < tolerance:
            break
            
    return reconstructedSignal, cosineComponents, nFunctionsUsed

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_rmsd(signal1: np.array, signal2: np.array) -> float:
    return np.sqrt(np.mean((signal1 - signal2) ** 2))
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def plot_array(x: np.array, y:np.array ):

    fig, ax = plt.subplots()
    ax.plot(x,y)

    plt.savefig("test.png")
    plt.close()



##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_signal_and_components(angle: np.array,
                                signal: np.array,
                                reconstructed_signal: np.array,
                                cosine_components: List[Tuple[float, np.array]],
                                out_dir: DirPath,
                                tag: str):
    ########
    plt.figure(figsize=(12, 6))

    angle = range(0,360, 10)
    ########

    plt.plot(angle, signal, label='QM Scan Energy', color='black',linewidth=2)
    ########
    for freq, cosine_component in cosine_components:
        plt.plot(angle, cosine_component, linestyle='--',
                 label=f'Cosine Component {freq:.2f}')
    ########
    plt.plot(angle, reconstructed_signal, label='Reconstructed Signal',
             color='red', linewidth=2, alpha=0.7)
    ########
    plt.xticks(np.arange(min(angle), max(angle) + 1, 60))
    plt.xlabel('Torsion Angle (Degrees)')
    plt.ylabel('Torsion Energy (kJ / mol)')
    plt.title(f'Fourier Tranform for {tag} Scan')
    plt.legend()
    plt.grid(True)
    ########
    save_png = p.join(out_dir, f"{tag}.png")
    plt.savefig(save_png)
    plt.close()
    ########


##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pad_energy_data(energyData: np.array, paddingFactor: int) -> np.array:
    return np.tile(energyData, paddingFactor)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pad_angle_data(angleData: np.array, paddingFactor: int) -> np.array:
    return np.tile(angleData, paddingFactor)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rescale_torsion_angles(angle):
    angle = angle % 360  # Reduce the angle to the 0-360 range
    if angle < 0:
        angle += 360  # Shift negative angles to the positive side
    return angle

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def perform_rfft(signal: np.array) -> np.array:
    return np.fft.rfft(signal)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_frequencies(signalLength: int, sampleSpacing: int) -> np.array:
    return np.fft.rfftfreq(signalLength, sampleSpacing)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def compute_amplitude_and_phase(fftResult: np.array, signalLength: int) -> Tuple[np.array, np.array]:
    amplitudes: np.array = np.abs(fftResult) * 2 / signalLength
    phases: np.array = np.angle(fftResult)

    return amplitudes, phases

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    
    outputDir = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/"
    torsionTopDir = p.join(outputDir, "02_torsion_scanning")
    torsionTag = 'O-C-CA-N'
    torsionDir = p.join(torsionTopDir, f"torsion_{torsionTag}")
    scanDir = p.join(torsionDir, "scan_data")
    energyCsv = p.join(scanDir, "final_scan_energies.csv")
    energyDf = pd.read_csv(energyCsv)

    fittingDir = p.join(torsionDir, "fitting")

    ##TODO: maxFunctions = 4
    
    fourier_transform_protocol(energyDf[torsionTag], torsionTag, fittingDir)
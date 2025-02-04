import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from os import path as p
from typing import Tuple, List

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
## dummy classes
class FilePath:
    pass
class DirPath:
    pass
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def dummy_inputs():
    outputDir = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/"
    torsionTag = 'O-C-CA-N'
    torsionDir = p.join(outputDir, f"torsion_{torsionTag}")
    scanDir = p.join(torsionDir, "scan_data")
    energyCsv = p.join(scanDir, "final_scan_energies.csv")
    energyDf = pd.read_csv(energyCsv)

    fittingDir = p.join(torsionDir, "fitting")

    return energyDf, torsionTag, fittingDir
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def main():
    energyDf, torsionTag, fittingDir = dummy_inputs()
    os.makedirs(fittingDir, exist_ok=True)

    energyData: np.array = energyDf[torsionTag].values

    energyDataPadded: np.array = pad_energy_data(energyData, paddingFactor=2)

    fit_scan_data(energyDataPadded, outDir=fittingDir, maxFunctions = 4)

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def perform_rfft(signal: np.array) -> np.array:
    return np.fft.rfft(signal)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_frequencies(signalLength: int, sampleSpacing: int) -> np.array:
    return np.fft.rfftfreq(signalLength, sampleSpacing) * 360
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def compute_amplitude_and_phase(fftResult: np.array, signalLength: int) -> Tuple[np.array, np.array]:
    amplitudes: np.array = np.abs(fftResult) * 2 / signalLength
    phases: np.array = np.angle(fftResult)
    ## round phases
    phases: np.array = np.degrees(phases)
    phases: np.array = np.round(phases).astype(int)
    phases: np.array = convert_phases_to_nearest(phases)
    return amplitudes, phases
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def convert_phases_to_nearest(phases: np.array) ->  np.array:
    return   np.where((180 - np.abs(phases)) < (np.abs(phases)), 180, 0)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def construct_cosine_components(dataDf: pd.DataFrame,
                                 angle: np.array,
                                     maxFunctions: int,
                                       tolerance: float = 1.0) -> Tuple[np.array, List[Tuple[float, np.array]], int]:
    # Initialize variables
    cosineComponents: list = []
    converged: bool = False
    previousSignal: np.array = np.zeros_like(angle)  # Initial previous signal

    while not converged:
        # Reset reconstructedSignal for each iteration
        reconstructedSignal: np.array = np.zeros_like(previousSignal)
        for index, row in dataDf.iterrows():
            print(index + 1)
            if index + 1 > maxFunctions:
                print("Max Functions reached!")
                nFunctionsUsed = maxFunctions
                converged = True
                break
            # Get variables from row
            A: float = row["Amplitude"]
            F: float = row["Frequency"]
            Phi: float = row["Phase"]
            # Construct cosine function using variables
            cosineComponent: np.array = A * np.cos(2 * np.pi * F * angle / 360 + np.radians(Phi))
            # Add to reconstructed signal
            reconstructedSignal += cosineComponent

            # Calculate RMSD between the current and previous reconstructed signals
            rmsd: float = calculate_rmsd(reconstructedSignal, previousSignal)

            # Check for convergence
            if rmsd < tolerance:
                print("reached Tol")
                previousSignal = reconstructedSignal.copy()
                cosineComponents.append((F, cosineComponent))
                nFunctionsUsed = index + 1
                converged = True
                break
            else:
                # Update previousSignal for the next iteration
                print("not reached tol")
                previousSignal = reconstructedSignal.copy()
                cosineComponents.append((F, cosineComponent))
    

    return previousSignal, cosineComponents, nFunctionsUsed
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_rmsd(signal1: np.array, signal2: np.array) -> float:
    return np.sqrt(np.mean((signal1 - signal2) ** 2))
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_signal_and_components(angle: np.array,
                                signal: np.array,
                                  reconstructedSignal: np.array,
                                    cosineComponents: List[Tuple[float, np.array]],
                                      outDir: DirPath,
                                        tag: str):
    ## set up a figure
    plt.figure(figsize=(12, 6))
    plt.plot(angle, signal, label='QM Scan Energy', color='black', linewidth=2)
    # plot the cosine components
    for freq, cosineComponent in cosineComponents:
        plt.plot(angle, cosineComponent, linestyle='--', label=f'Cosine Component {freq:.2f}')
    # plot reconstructed signal
    plt.plot(angle, reconstructedSignal, label='Reconstructed Signal', color='red', linewidth=2, alpha=0.7)
    # add labels and title
    plt.xlabel('angle')
    plt.ylabel('Amplitude')
    plt.title('Original Signal and Cosine Components')
    plt.legend()
    plt.grid(True)
    # save the figure and close
    savePng = p.join(outDir, f"{tag}.png")
    plt.savefig(savePng)
    plt.close()
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def fit_scan_data(signal: np.array,
                   outDir: DirPath,
                     sampleSpacing: float = 10.0,
                       tag: str = 'cosine_components',
                         maxFunctions: int = 10) -> None:
    
    ## calculate signal length
    signalLength: int = len(signal)
    ## run reverse fourier transform
    fftResult: np.array = perform_rfft(signal)
    ## get frequencies, amplitudes and phases
    frequencies: np.array = get_frequencies(signalLength, sampleSpacing)
    amplitudes, phases = compute_amplitude_and_phase(fftResult, signalLength)
    ## construct angle x-axis
    angle: np.array = np.arange(signalLength) * sampleSpacing
    ## convert data to dataframe
    dataDf: pd.DataFrame = convert_data_to_df(frequencies, amplitudes, phases)
    ## construct cosine components from parameters
    reconstructedSignal, cosineComponents, nFunctionsUsed = construct_cosine_components(dataDf, angle, maxFunctions)
    ## plot signal and components
    plot_signal_and_components(angle, signal, reconstructedSignal, cosineComponents, outDir, tag)
    ## write data to csv file
    outCsv: FilePath = p.join(outDir, f"{tag}.csv")
    dataDf.iloc[:nFunctionsUsed].to_csv(outCsv)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def convert_data_to_df(frequencies: np.array, amplitudes: np.array, phases: np.array) -> pd.DataFrame:
    data = {"Frequency": frequencies,
            "Amplitude": amplitudes,
            "Phase": phases}    
    ## convert to df
    dataDf: pd.DataFrame = pd.DataFrame(data)

    # remove first component (this will be a constant)
    dataDf = dataDf.iloc[1:]

    # sort data by amplitude
    dataDf = dataDf.sort_values(by='Amplitude', ascending=False).reset_index(drop=True)

    return dataDf
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pad_energy_data(energyData: np.array, paddingFactor: int) -> np.array:
    energyDataPadded = np.concatenate([energyData[::-1] if i % 2 else energyData for i in range(paddingFactor * 2)])
    return energyDataPadded
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pad_angle_data(angleData: np.array, paddingFactor: int) -> np.array:
    angleDataPadded = np.concatenate([angleData + 360 * i for i in range(paddingFactor * 2)])
    return angleDataPadded
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    main()
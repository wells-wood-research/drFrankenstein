import numpy as np
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
def dummy_inputs():
    outputDir = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/"
    torsionTopDir = p.join(outputDir, "02_torsion_scanning")
    torsionTag = 'O-C-CA-N'
    torsionDir = p.join(torsionTopDir, f"torsion_{torsionTag}")
    scanDir = p.join(torsionDir, "scan_data")
    energyCsv = p.join(scanDir, "final_scan_energies.csv")
    energyDf = pd.read_csv(energyCsv)

    fittingDir = p.join(torsionDir, "fitting")

    ##TODO: maxFunctions = 4

    return energyDf, torsionTag, fittingDir
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def main():
    energyDf, torsionTag, fittingDir = dummy_inputs()
    os.makedirs(fittingDir, exist_ok=True)
    ## load energy data [range of -170 to 180 degrees]
    ## convert to a [0 - 360] degree range
    energyDf["Angle"] = energyDf["Angle"].apply(rescale_torsion_angles)
    energyDf = energyDf.sort_values(by="Angle", ascending=True)

    energyData: np.array = energyDf[torsionTag].values

    energyDataPadded: np.array = pad_energy_data(energyData, paddingFactor=1)

    fit_scan_data(energyDataPadded, outDir=fittingDir, maxFunctions = 4, tag = torsionTag)

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
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
    return np.fft.rfftfreq(signalLength, sampleSpacing) * 360
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def compute_amplitude_and_phase(fftResult: np.array, signalLength: int) -> Tuple[np.array, np.array]:
    amplitudes: np.array = np.abs(fftResult) * 2 / signalLength
    phases: np.array = np.angle(fftResult)
    ## round phases
    phases: np.array = np.degrees(phases)
    phases: np.array = np.round(phases).astype(int)
    # phases: np.array = convert_phases_to_nearest(phases)
    return amplitudes, phases
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def convert_phases_to_nearest(phases: np.array) ->  np.array:
    return   np.where((180 - np.abs(phases)) < (np.abs(phases)), 180, 0)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def construct_cosine_components(dataDf: pd.DataFrame,
                                 angle: np.array,
                                     maxFunctions: int,
                                       tolerance: float = 0.5) -> Tuple[np.array, List[Tuple[float, np.array]], int]:
    # Initialize variables
    cosineComponents: list = []
    converged: bool = False
    previousSignal: np.array = np.zeros_like(angle)  # Initial previous signal

    while not converged:
        # Reset reconstructedSignal for each iteration
        reconstructedSignal: np.array = np.zeros_like(previousSignal)
        for index, row in dataDf.iterrows():
            print(index + 1)
            if index + 1  == maxFunctions:
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
def plot_signal_and_components_html(angle: np.array,
                               signal: np.array,
                               cosineComponents: List[Tuple[float, np.array]],
                               outDir: DirPath,
                               tag: str):
    fig = go.Figure()

    fig.add_trace(go.Scatter(x=angle, y=signal, mode='lines',
                             name='QM Scan Energy', line=dict(color='black', width=2)))

    # Initialize reconstructed signal as zeros
    reconstructedSignal = np.zeros_like(angle)

    for freq, cosineComponent in cosineComponents:
        fig.add_trace(go.Scatter(x=angle, y=cosineComponent, mode='lines',
                                 name=f'Cosine Component {freq:.2f}', 
                                 line=dict(dash='dash')))
        # Add each cosine component to the reconstructed signal
        reconstructedSignal += cosineComponent

    fig.add_trace(go.Scatter(x=angle, y=reconstructedSignal, mode='lines',
                             name='Reconstructed Signal', 
                             line=dict(color='red', width=2), opacity=0.7))

    fig.update_layout(
        xaxis_title='angle',
        yaxis_title='Amplitude',
        title='Original Signal and Cosine Components',
        legend=dict(x=0, y=1),
        xaxis=dict(tickmode='array', tickvals=np.arange(min(angle), max(angle) + 1, 60)),
        template='plotly_white'
    )

    saveHtml = p.join(outDir, f"{tag}.html")
    fig.write_html(saveHtml)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_signal_and_components(angle: np.array,
                                signal: np.array,
                                reconstructed_signal: np.array,
                                cosine_components: List[Tuple[float, np.array]],
                                out_dir: DirPath,
                                tag: str):
    ########
    plt.figure(figsize=(12, 6))
    angle_mask = (angle >= 0) & (angle <= 360)
    angle = angle[angle_mask]
    signal = signal[angle_mask]
    reconstructed_signal = reconstructed_signal[angle_mask]
    ########
    plt.plot(angle, signal, label='QM Scan Energy', color='black',
             linewidth=2)
    ########
    for freq, cosine_component in cosine_components:
        cosine_component = cosine_component[angle_mask]
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
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def fit_scan_data(signal: np.array,
                   outDir: DirPath,
                     maxFunctions: int,
                     sampleSpacing: float = 10.0,
                       tag: str = 'cosine_components') -> None:
    
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
    print(dataDf)
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

    # # remove first component (this will be a constant)
    dataDf = dataDf.iloc[1:]

    # sort data by amplitude
    dataDf = dataDf.sort_values(by='Amplitude', ascending=False).reset_index(drop=True)

    return dataDf
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pad_energy_data(energyData: np.array, paddingFactor: int) -> np.array:
    energyDataPadded = np.concatenate([energyData[::-1] if i % 2 else energyData for i in range(paddingFactor * 2)])

    energyDataPadded = np.concatenate([energyData for I in range(paddingFactor * 2)])

    return energyDataPadded
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pad_angle_data(angleData: np.array, paddingFactor: int) -> np.array:
    angleDataPadded = np.concatenate([angleData + 360 * i for i in range(paddingFactor * 2)])
    return angleDataPadded
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    main()
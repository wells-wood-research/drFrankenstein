import numpy as np

import matplotlib
matplotlib.use('Agg')  # Set the Agg backend before importing pyplot
import pandas as pd
from os import path as p
from typing import Tuple, List

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
## dummy classes
class FilePath:
    pass
class DirPath:
    pass

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def fourier_transform_protocol(qmTorsionEnergy, torsionTag, torsionFittingDir, sampleSpacing=10, maxFunctions=3, forceField = "AMBER", l2Damping = 0.1):    
    energyDataPadded: np.array = pad_energy_data(qmTorsionEnergy, paddingFactor=3)
    ## calculate signal length
    signalLength: int = len(energyDataPadded)
    ## run reverse fourier transform
    fftResult: np.array = perform_rfft(energyDataPadded)
    ## get frequencies, amplitudes and phases
    frequencies: np.array = get_frequencies(signalLength, sampleSpacing)
    amplitudes, phases = compute_amplitude_and_phase(fftResult, signalLength)

    amplitudes = apply_l2_damping(amplitudes, l2Damping)

    ## construct angle x-axis
    angle: np.array = np.arange(signalLength) * sampleSpacing
    ## convert data to dataframe
    fourierDf: pd.DataFrame = convert_fourier_params_to_df(frequencies, amplitudes, phases)
    paramDf: pd.DataFrame = convert_params_to_amber_charmm_format(fourierDf)
    if forceField == "AMBER":
        ## construct cosine components from parameters
        reconstructedSignal, cosineComponents, nFunctionsUsed = construct_cosine_components_AMBER(paramDf, angle, maxFunctions)

    elif forceField == "CHARMM":
        reconstructedSignal, cosineComponents, nFunctionsUsed = construct_cosine_components_CHARMM(paramDf, angle, maxFunctions)

        ## write data to csv file
    outCsv: FilePath = p.join(torsionFittingDir, f"{torsionTag}.csv")
    paramDf.iloc[:nFunctionsUsed].to_csv(outCsv)

    return paramDf.iloc[:nFunctionsUsed], cosineComponents   
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def apply_l2_damping(amplitudes: np.array, l2Damping: float) -> np.array:
    """
    Applies a L2 damping to the amplitudes
    This may help prevent escalating amplitudes as nShuffles increases

    Args:
        amplitudes (np.array): array of amplitudes
        l2Damping (float): damping factor

    Returns:
        dampenedAmplitudes (np.array): dampened amplitudes
    
    """
    dampenedAmplitudes = amplitudes / (1 + l2Damping * np.abs(amplitudes))

    return dampenedAmplitudes


def convert_fourier_params_to_df(frequencies: np.array, amplitudes: np.array, phases: np.array) -> pd.DataFrame:
    data = {"Frequency": frequencies, "Amplitude": amplitudes, "Phase": phases}
    dataDf = pd.DataFrame(data)
    dataDf = dataDf.sort_values(by='Amplitude', ascending=False).reset_index(drop=True)
    return dataDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def convert_params_to_amber_charmm_format(fourierDf: pd.DataFrame) -> pd.DataFrame:
    paramDf = pd.DataFrame()
    paramDf["Amplitude"] = fourierDf["Amplitude"]
    paramDf["Period"] = np.degrees(2 * np.pi * fourierDf["Frequency"])
    paramDf["Phase"] = np.degrees(fourierDf["Phase"]) * -1

    ## remove period == 0 (DC) Signal
    paramDf = paramDf[paramDf["Period"] > 0]

    paramDf.sort_values(by="Amplitude", ascending=False, inplace=True,ignore_index=True)
    return paramDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def convert_params_to_charmm_format(fourierDf: pd.DataFrame) -> pd.DataFrame:
    charmmDf = pd.DataFrame()
    charmmDf["Amplitude"] = fourierDf["Amplitude"]
    charmmDf["Period"] = np.degrees(2 * np.pi * fourierDf["Frequency"])
    charmmDf["Phase"] = np.degrees(fourierDf["Phase"]) * -1
    ## remove period == 0 (DC) Signal
    charmmDf = charmmDf[charmmDf["Period"] > 0]
    charmmDf.sort_values(by="Amplitude", ascending=False, inplace=True, ignore_index=True)

    return charmmDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def construct_cosine_components_CHARMM(charmmParamDf: pd.DataFrame,
                                 angle: np.array,
                                     maxFunctions: int,
                                       tolerance: float = 0.2) -> Tuple[np.array, List[Tuple[float, np.array]], int]:
    
    amplitudes = charmmParamDf["Amplitude"]
    periods = charmmParamDf["Period"]
    phases = charmmParamDf["Phase"]

    sampleSpacing = 10
    signalLength = 36
    
    # Angle array
    angle = np.arange(signalLength) * sampleSpacing  # Shape: (N,)
    
    # Initialize reconstructed signal
    reconstructedSignal = np.zeros(signalLength)
    previousSignal = np.zeros(signalLength)

    ## init mean average error
    meanAverageError = np.inf
    ## collect cosine components for plotting
    cosineComponents = []
    # Construct each cosine component
    while True:
        for nFunctionsUsed in range(1, maxFunctions + 1):
            charmmComponent = amplitudes[nFunctionsUsed] * (1 + np.cos(np.radians(periods[nFunctionsUsed] * angle - phases[nFunctionsUsed])))
            previousSignal = reconstructedSignal.copy()
            reconstructedSignal += charmmComponent
            meanAverageError = np.mean(np.abs(reconstructedSignal - previousSignal))
            cosineComponents.append((nFunctionsUsed, charmmComponent))
            if meanAverageError < tolerance:
                break
        break  # Exit the while loop after the for loop completes or breaks

    return reconstructedSignal, cosineComponents, nFunctionsUsed
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def construct_cosine_components_AMBER(amberParamDf: pd.DataFrame,
                                 angle: np.array,
                                     maxFunctions: int,
                                       tolerance: float = 0.05) -> Tuple[np.array, List[Tuple[float, np.array]], int]:
    

    amplitudes = amberParamDf["Amplitude"]
    phases = amberParamDf["Phase"]
    periods = amberParamDf["Period"]


    sampleSpacing = 10
    signalLength = 36
    
    # Angle array
    angle = np.arange(signalLength) * sampleSpacing  # Shape: (N,)
    
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
        # if meanAverageError < tolerance:
        break
    return reconstructedSignal, cosineComponents, nFunctionsUsed

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_rmsd(signal1: np.array, signal2: np.array) -> float:
    return np.sqrt(np.mean((signal1 - signal2) ** 2))
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
    
    raise NotImplementedError
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
from . import Stitching_Assistant


def fourier_transform_protocol(qmTorsionEnergy, torsionTag, torsionFittingDir, config):
                            #    sampleSpacing=10, maxFunctions=3, forceField = "AMBER", l2Damping = 0.1):    

    ## unpack config
    maxCosineFunctions = config["runtimeInfo"]["madeByStitching"]["maxTorsions"]
    l2Damping = config["parameterFittingInfo"]["l2DampingFactor"]
    forceField = config["parameterFittingInfo"]["forceField"]
    sampleSpacing = 10
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
        selectedParamDf, cosineComponents, nFunctionsUsed = construct_cosine_components_AMBER(
            paramDf, angle, maxCosineFunctions, qmTorsionEnergy
        )

    elif forceField == "CHARMM":
        selectedParamDf, cosineComponents, nFunctionsUsed = construct_cosine_components_CHARMM(
            paramDf, angle, maxCosineFunctions, qmTorsionEnergy
        )

    ## write data to csv file
    outCsv: FilePath = p.join(torsionFittingDir, f"{torsionTag}.csv")
    selectedParamDf.to_csv(outCsv)

    reconstructedSignal = np.sum([component for _, component in cosineComponents], axis=0) if cosineComponents else np.zeros(len(qmTorsionEnergy))
    ## calculate mean average error
    meanAverageError = np.mean(np.abs(reconstructedSignal - qmTorsionEnergy))
    fitScore = Stitching_Assistant.calculate_profile_fit_score(qmTorsionEnergy, reconstructedSignal)

    return selectedParamDf, cosineComponents, reconstructedSignal, meanAverageError, fitScore   
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
def construct_cosine_components_CHARMM(
    charmmParamDf: pd.DataFrame,
    angle: np.array,
    maxFunctions: int,
    qmTorsionEnergy: np.array,
    tolerance: float = 0.2,
) -> Tuple[np.array, List[Tuple[float, np.array]], int]:
    return _construct_cosine_components(charmmParamDf, maxFunctions, qmTorsionEnergy)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def construct_cosine_components_AMBER(
    amberParamDf: pd.DataFrame,
    angle: np.array,
    maxFunctions: int,
    qmTorsionEnergy: np.array,
    tolerance: float = 0.05,
) -> Tuple[np.array, List[Tuple[float, np.array]], int]:
    return _construct_cosine_components(amberParamDf, maxFunctions, qmTorsionEnergy)


def _construct_cosine_components(
    paramDf: pd.DataFrame,
    maxFunctions: int,
    qmTorsionEnergy: np.array,
) -> Tuple[pd.DataFrame, List[Tuple[float, np.array]], int]:
    sampleSpacing = 10
    signalLength = len(qmTorsionEnergy)
    angle = np.arange(signalLength) * sampleSpacing

    reconstructedSignal = np.zeros(signalLength)
    bestScore = Stitching_Assistant.calculate_profile_fit_score(qmTorsionEnergy, reconstructedSignal)
    cosineComponents = []
    selectedRowIndexes = []

    for _, row in paramDf.iloc[:maxFunctions].iterrows():
        candidateComponent = row["Amplitude"] * (1 + np.cos(np.radians(row["Period"] * angle - row["Phase"])))
        candidateSignal = reconstructedSignal + candidateComponent
        candidateScore = Stitching_Assistant.calculate_profile_fit_score(qmTorsionEnergy, candidateSignal)
        if candidateScore < bestScore:
            reconstructedSignal = candidateSignal
            bestScore = candidateScore
            selectedRowIndexes.append(row.name)
            cosineComponents.append((len(cosineComponents) + 1, candidateComponent))

    selectedParamDf = paramDf.loc[selectedRowIndexes]
    return selectedParamDf, cosineComponents, len(cosineComponents)

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

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
def fourier_transform_protocol(qm_torsion_energy: np.ndarray, 
                                torsion_tag: str, 
                                torsion_fitting_dir: DirPath, 
                                sample_spacing: int =10, 
                                max_functions: int =3, 
                                force_field: str = "AMBER", 
                                l2_damping: float = 0.1) -> Tuple[pd.DataFrame, dict]:
    """
    Performs Fourier transform on QM torsion energy to derive force field parameters.

    Args:
        qm_torsion_energy: NumPy array of QM torsion energies.
        torsion_tag: Identifier for the torsion.
        torsion_fitting_dir: Directory to save output CSV.
        sample_spacing: Spacing between samples in the energy profile.
        max_functions: Maximum number of cosine functions to use for reconstruction.
        force_field: The target force field ('AMBER' or 'CHARMM').
        l2_damping: L2 damping factor for amplitudes.

    Returns:
        A tuple containing a DataFrame of the derived parameters (up to nFunctionsUsed)
        and a dictionary of the cosine components used in the reconstruction.
    """
    energyDataPadded: np.ndarray = pad_energy_data(energy_data=qm_torsion_energy, padding_factor=3)
    ## calculate signal length
    signalLength: int = len(energyDataPadded)
    ## run reverse fourier transform
    fftResult: np.ndarray = perform_rfft(signal=energyDataPadded)
    ## get frequencies, amplitudes and phases
    frequencies: np.ndarray = get_frequencies(signal_length=signalLength, sample_spacing=sample_spacing)
    amplitudes, phases = compute_amplitude_and_phase(fft_result=fftResult, signal_length=signalLength)

    amplitudes = apply_l2_damping(amplitudes=amplitudes, l2_damping=l2_damping)

    ## construct angle x-axis
    angle: np.ndarray = np.arange(signalLength) * sample_spacing
    ## convert data to dataframe
    fourierDf: pd.DataFrame = convert_fourier_params_to_df(frequencies=frequencies, amplitudes=amplitudes, phases=phases)
    paramDf: pd.DataFrame = convert_params_to_amber_charmm_format(fourier_df=fourierDf)
    if force_field == "AMBER":
        ## construct cosine components from parameters
        reconstructedSignal, cosineComponents, nFunctionsUsed = construct_cosine_components_AMBER(amber_param_df=paramDf, angle=angle, max_functions=max_functions)

    elif force_field == "CHARMM":
        reconstructedSignal, cosineComponents, nFunctionsUsed = construct_cosine_components_CHARMM(charmm_param_df=paramDf, angle=angle, max_functions=max_functions)
    else:
        raise ValueError(f"Unsupported force_field: {force_field}. Must be 'AMBER' or 'CHARMM'.")

        ## write data to csv file
    outCsv: FilePath = p.join(torsion_fitting_dir, f"{torsion_tag}.csv") # type: ignore
    paramDf.iloc[:nFunctionsUsed].to_csv(outCsv)

    return paramDf.iloc[:nFunctionsUsed], cosineComponents   
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def apply_l2_damping(amplitudes: np.ndarray, l2_damping: float) -> np.ndarray:
    """
    Applies a L2 damping to the amplitudes
    This may help prevent escalating amplitudes as nShuffles increases

    Args:
        amplitudes (np.array): array of amplitudes
        l2Damping (float): damping factor

    Returns:
        dampenedAmplitudes (np.array): dampened amplitudes
    
    """
    dampenedAmplitudes = amplitudes / (1 + l2_damping * np.abs(amplitudes))

    return dampenedAmplitudes


def convert_fourier_params_to_df(frequencies: np.ndarray, amplitudes: np.ndarray, phases: np.ndarray) -> pd.DataFrame:
    data = {"Frequency": frequencies, "Amplitude": amplitudes, "Phase": phases}
    dataDf = pd.DataFrame(data)
    dataDf = dataDf.sort_values(by='Amplitude', ascending=False).reset_index(drop=True)
    return dataDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def convert_params_to_amber_charmm_format(fourier_df: pd.DataFrame) -> pd.DataFrame:
    paramDf = pd.DataFrame()
    paramDf["Amplitude"] = fourier_df["Amplitude"]
    paramDf["Period"] = np.degrees(2 * np.pi * fourier_df["Frequency"])
    paramDf["Phase"] = np.degrees(fourier_df["Phase"]) * -1

    ## remove period == 0 (DC) Signal
    paramDf = paramDf[paramDf["Period"] > 0]

    paramDf.sort_values(by="Amplitude", ascending=False, inplace=True,ignore_index=True)
    return paramDf
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def convert_params_to_charmm_format(fourier_df: pd.DataFrame) -> pd.DataFrame:
    charmmDf = pd.DataFrame()
    charmmDf["Amplitude"] = fourier_df["Amplitude"]
    charmmDf["Period"] = np.degrees(2 * np.pi * fourier_df["Frequency"])
    charmmDf["Phase"] = np.degrees(fourier_df["Phase"]) * -1
    ## remove period == 0 (DC) Signal
    charmmDf = charmmDf[charmmDf["Period"] > 0]
    charmmDf.sort_values(by="Amplitude", ascending=False, inplace=True, ignore_index=True)

    return charmmDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def construct_cosine_components_CHARMM(charmm_param_df: pd.DataFrame,
                                 angle: np.ndarray,
                                     max_functions: int,
                                       tolerance: float = 0.2) -> Tuple[np.ndarray, List[Tuple[float, np.ndarray]], int]:
    
    amplitudes = charmm_param_df["Amplitude"]
    periods = charmm_param_df["Period"] # Corrected: Access 'Period' from DataFrame
    phases = charmm_param_df["Phase"]

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
        for nFunctionsUsed_idx in range(max_functions): # Iterate up to max_functions -1 (for 0-indexed df)
            # Ensure index is within bounds of the DataFrame
            if nFunctionsUsed_idx >= len(amplitudes):
                break # Not enough functions in DataFrame as requested by max_functions

            amplitude_val = amplitudes.iloc[nFunctionsUsed_idx]
            period_val = periods.iloc[nFunctionsUsed_idx]
            phase_val = phases.iloc[nFunctionsUsed_idx]

            charmmComponent = amplitude_val * (1 + np.cos(np.radians(period_val * angle - phase_val)))
            previousSignal = reconstructedSignal.copy()
            reconstructedSignal += charmmComponent
            meanAverageError = np.mean(np.abs(reconstructedSignal - previousSignal))
            # Store with 1-based indexing for nFunctionsUsed if that's the intended meaning
            cosineComponents.append((nFunctionsUsed_idx + 1, charmmComponent))
            if meanAverageError < tolerance and nFunctionsUsed_idx > 0: # Avoid breaking on first component if MAE is low
                break
        nFunctionsUsed = nFunctionsUsed_idx + 1 # Actual number of functions used
        break  # Exit the while loop after the for loop completes or breaks

    return reconstructedSignal, cosineComponents, nFunctionsUsed
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def construct_cosine_components_AMBER(amber_param_df: pd.DataFrame,
                                 angle: np.ndarray,
                                     max_functions: int,
                                       tolerance: float = 0.05) -> Tuple[np.ndarray, List[Tuple[float, np.ndarray]], int]:
    
    amplitudes = amber_param_df["Amplitude"]
    periods = amber_param_df["Period"] # Corrected: Access 'Period' from DataFrame
    phases = amber_param_df["Phase"]

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
        for nFunctionsUsed_idx in range(max_functions): # Iterate up to max_functions -1 (for 0-indexed df)
             # Ensure index is within bounds of the DataFrame
            if nFunctionsUsed_idx >= len(amplitudes):
                break # Not enough functions in DataFrame as requested by max_functions

            amplitude_val = amplitudes.iloc[nFunctionsUsed_idx]
            period_val = periods.iloc[nFunctionsUsed_idx]
            phase_val = phases.iloc[nFunctionsUsed_idx]
            
            amberComponent =  amplitude_val * (1 + np.cos(np.radians(period_val * angle - phase_val)))

            previousSignal = reconstructedSignal.copy()
            reconstructedSignal += amberComponent
            meanAverageError =  np.mean(np.abs(reconstructedSignal - previousSignal))
            # Store with 1-based indexing for nFunctionsUsed if that's the intended meaning
            cosineComponents.append((nFunctionsUsed_idx + 1 , amberComponent))
            if meanAverageError < tolerance and nFunctionsUsed_idx > 0: # Avoid breaking on first component
                break
        nFunctionsUsed = nFunctionsUsed_idx + 1 # Actual number of functions used
        break
    return reconstructedSignal, cosineComponents, nFunctionsUsed

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_rmsd(signal1: np.ndarray, signal2: np.ndarray) -> float:
    return np.sqrt(np.mean((signal1 - signal2) ** 2))
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def pad_energy_data(energy_data: np.ndarray, padding_factor: int) -> np.ndarray:
    return np.tile(energy_data, padding_factor)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pad_angle_data(angle_data: np.ndarray, padding_factor: int) -> np.ndarray:
    return np.tile(angle_data, padding_factor)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rescale_torsion_angles(angle: float) -> float: # Assuming float for single angle, could be np.ndarray
    angle = angle % 360  # Reduce the angle to the 0-360 range
    if angle < 0:
        angle += 360  # Shift negative angles to the positive side
    return angle

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def perform_rfft(signal: np.ndarray) -> np.ndarray:
    return np.fft.rfft(signal)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_frequencies(signal_length: int, sample_spacing: int) -> np.ndarray:
    return np.fft.rfftfreq(signal_length, sample_spacing)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def compute_amplitude_and_phase(fft_result: np.ndarray, signal_length: int) -> Tuple[np.ndarray, np.ndarray]:
    amplitudes: np.ndarray = np.abs(fft_result) * 2 / signal_length
    phases: np.ndarray = np.angle(fft_result)

    return amplitudes, phases

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    
    raise NotImplementedError
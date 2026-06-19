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


def _coerce_positive_float(value: object, default: float) -> float:
    """Parse a positive float with fallback."""
    try:
        parsed = float(value)
        if parsed <= 0:
            return default
        return parsed
    except Exception:
        return default


def fourier_transform_protocol(
    qmTorsionEnergy: np.array,
    torsionTag: str,
    torsionFittingDir: DirPath,
    config: dict,
    qmTotalEnergy: np.ndarray | None = None,
    mmTotalEnergy: np.ndarray | None = None,
    mmTorsionEnergy: np.ndarray | None = None,
):
                            #    sampleSpacing=10, maxFunctions=3, forceField = "AMBER", l2Damping = 0.1):
    """Fit Fourier torsion parameters and return the selected components."""

    ## unpack config
    maxCosineFunctions = config["runtimeInfo"]["madeByStitching"]["maxTorsions"]
    l2Damping = config["parameterFittingInfo"]["l2DampingFactor"]
    forceField = config["parameterFittingInfo"]["forceField"]
    acceptableFitScore = config["parameterFittingInfo"].get("converganceTolerance", None)
    cosineComplexityPenalty = _coerce_positive_float(
        config["parameterFittingInfo"].get("cosineComplexityPenalty", 0.03), 0.03
    )
    minScoreImprovement = _coerce_positive_float(
        config["parameterFittingInfo"].get("cosineMinScoreImprovement", 1e-4), 1e-4
    )
    hasTotalObjective = (
        qmTotalEnergy is not None
        and mmTotalEnergy is not None
        and mmTorsionEnergy is not None
        and len(qmTotalEnergy) == len(qmTorsionEnergy) == len(mmTotalEnergy) == len(mmTorsionEnergy)
    )
    if hasTotalObjective:
        torsionObjectiveWeight = _coerce_positive_float(
            config["parameterFittingInfo"].get("torsionObjectiveWeight", 0.5), 0.5
        )
        totalObjectiveWeight = _coerce_positive_float(
            config["parameterFittingInfo"].get("totalObjectiveWeight", 0.5), 0.5
        )
    else:
        torsionObjectiveWeight = 1.0
        totalObjectiveWeight = 0.0
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
            paramDf,
            angle,
            maxCosineFunctions,
            qmTorsionEnergy,
            qmTotalEnergy=qmTotalEnergy,
            mmTotalEnergy=mmTotalEnergy,
            mmTorsionEnergy=mmTorsionEnergy,
            torsionObjectiveWeight=torsionObjectiveWeight,
            totalObjectiveWeight=totalObjectiveWeight,
            acceptableFitScore=acceptableFitScore,
            complexityPenalty=cosineComplexityPenalty,
            minScoreImprovement=minScoreImprovement,
        )

    elif forceField == "CHARMM":
        selectedParamDf, cosineComponents, nFunctionsUsed = construct_cosine_components_CHARMM(
            paramDf,
            angle,
            maxCosineFunctions,
            qmTorsionEnergy,
            qmTotalEnergy=qmTotalEnergy,
            mmTotalEnergy=mmTotalEnergy,
            mmTorsionEnergy=mmTorsionEnergy,
            torsionObjectiveWeight=torsionObjectiveWeight,
            totalObjectiveWeight=totalObjectiveWeight,
            acceptableFitScore=acceptableFitScore,
            complexityPenalty=cosineComplexityPenalty,
            minScoreImprovement=minScoreImprovement,
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
    """Apply L2 damping to Fourier amplitudes."""
    """
    Applies L2 (ridge-style) damping to amplitudes.

    Scales amplitudes by 1 / (1 + l2Damping), independent of amplitude magnitude.

    Args:
        amplitudes (np.array): array of amplitudes
        l2Damping (float): damping factor

    Returns:
        dampenedAmplitudes (np.array): dampened amplitudes
    """
    # If no damping requested, return original amplitudes unchanged
    try:
        ld = float(l2Damping)
    except Exception:
        ld = 0.0
    if ld == 0.0:
        return amplitudes

    dampenedAmplitudes = amplitudes / (1.0 + ld)

    return dampenedAmplitudes


def convert_fourier_params_to_df(frequencies: np.array, amplitudes: np.array, phases: np.array) -> pd.DataFrame:
    """Combine Fourier arrays into a sorted DataFrame."""
    data = {"Frequency": frequencies, "Amplitude": amplitudes, "Phase": phases}
    dataDf = pd.DataFrame(data)
    dataDf = dataDf.sort_values(by='Amplitude', ascending=False).reset_index(drop=True)
    return dataDf

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def convert_params_to_amber_charmm_format(fourierDf: pd.DataFrame) -> pd.DataFrame:
    """Convert Fourier parameters into AMBER/CHARMM torsion format."""
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
    """Convert Fourier parameters into CHARMM torsion format."""
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
    qmTotalEnergy: np.ndarray | None = None,
    mmTotalEnergy: np.ndarray | None = None,
    mmTorsionEnergy: np.ndarray | None = None,
    torsionObjectiveWeight: float = 1.0,
    totalObjectiveWeight: float = 0.0,
    acceptableFitScore: float | None = None,
    complexityPenalty: float = 0.03,
    minScoreImprovement: float = 1e-4,
    tolerance: float = 0.2,
) -> Tuple[np.array, List[Tuple[float, np.array]], int]:
    """Construct cosine components for CHARMM torsions."""
    return _construct_cosine_components(
        charmmParamDf,
        maxFunctions,
        qmTorsionEnergy,
        qmTotalEnergy=qmTotalEnergy,
        mmTotalEnergy=mmTotalEnergy,
        mmTorsionEnergy=mmTorsionEnergy,
        torsionObjectiveWeight=torsionObjectiveWeight,
        totalObjectiveWeight=totalObjectiveWeight,
        acceptableFitScore=acceptableFitScore,
        complexityPenalty=complexityPenalty,
        minScoreImprovement=minScoreImprovement,
    )
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def construct_cosine_components_AMBER(
    amberParamDf: pd.DataFrame,
    angle: np.array,
    maxFunctions: int,
    qmTorsionEnergy: np.array,
    qmTotalEnergy: np.ndarray | None = None,
    mmTotalEnergy: np.ndarray | None = None,
    mmTorsionEnergy: np.ndarray | None = None,
    torsionObjectiveWeight: float = 1.0,
    totalObjectiveWeight: float = 0.0,
    acceptableFitScore: float | None = None,
    complexityPenalty: float = 0.03,
    minScoreImprovement: float = 1e-4,
    tolerance: float = 0.05,
) -> Tuple[np.array, List[Tuple[float, np.array]], int]:
    """Construct cosine components for AMBER torsions."""
    return _construct_cosine_components(
        amberParamDf,
        maxFunctions,
        qmTorsionEnergy,
        qmTotalEnergy=qmTotalEnergy,
        mmTotalEnergy=mmTotalEnergy,
        mmTorsionEnergy=mmTorsionEnergy,
        torsionObjectiveWeight=torsionObjectiveWeight,
        totalObjectiveWeight=totalObjectiveWeight,
        acceptableFitScore=acceptableFitScore,
        complexityPenalty=complexityPenalty,
        minScoreImprovement=minScoreImprovement,
    )


def _score_candidate(
    candidateTorsionSignal: np.ndarray,
    qmTorsionEnergy: np.ndarray,
    qmTotalEnergy: np.ndarray | None = None,
    mmTotalEnergy: np.ndarray | None = None,
    mmTorsionEnergy: np.ndarray | None = None,
    torsionObjectiveWeight: float = 1.0,
    totalObjectiveWeight: float = 0.0,
) -> tuple[float, float | None, float]:
    """Return (torsion score, total score or None, blended score)."""
    torsionScore = Stitching_Assistant.calculate_profile_fit_score(qmTorsionEnergy, candidateTorsionSignal)
    hasTotalObjective = (
        qmTotalEnergy is not None
        and mmTotalEnergy is not None
        and mmTorsionEnergy is not None
        and len(qmTotalEnergy) == len(candidateTorsionSignal) == len(mmTotalEnergy) == len(mmTorsionEnergy)
    )

    if not hasTotalObjective:
        return torsionScore, None, torsionScore

    candidateTotalSignal = mmTotalEnergy - mmTorsionEnergy + candidateTorsionSignal
    totalScore = Stitching_Assistant.calculate_profile_fit_score(qmTotalEnergy, candidateTotalSignal)
    weightSum = torsionObjectiveWeight + totalObjectiveWeight
    if weightSum <= 0:
        weightSum = 1.0
        torsionObjectiveWeight = 1.0
        totalObjectiveWeight = 0.0
    blendedScore = (
        torsionObjectiveWeight * torsionScore + totalObjectiveWeight * totalScore
    ) / weightSum
    return torsionScore, totalScore, blendedScore


def _construct_cosine_components(
    paramDf: pd.DataFrame,
    maxFunctions: int,
    qmTorsionEnergy: np.array,
    qmTotalEnergy: np.ndarray | None = None,
    mmTotalEnergy: np.ndarray | None = None,
    mmTorsionEnergy: np.ndarray | None = None,
    torsionObjectiveWeight: float = 1.0,
    totalObjectiveWeight: float = 0.0,
    acceptableFitScore: float | None = None,
    complexityPenalty: float = 0.03,
    minScoreImprovement: float = 1e-4,
) -> Tuple[pd.DataFrame, List[Tuple[float, np.array]], int]:
    """Select a compact cosine set that reaches acceptable fit when possible."""
    sampleSpacing = 10
    signalLength = len(qmTorsionEnergy)
    angle = np.arange(signalLength) * sampleSpacing

    reconstructedSignal = np.zeros(signalLength)
    bestTorsionScore, bestTotalScore, bestBlendedScore = _score_candidate(
        reconstructedSignal,
        qmTorsionEnergy,
        qmTotalEnergy=qmTotalEnergy,
        mmTotalEnergy=mmTotalEnergy,
        mmTorsionEnergy=mmTorsionEnergy,
        torsionObjectiveWeight=torsionObjectiveWeight,
        totalObjectiveWeight=totalObjectiveWeight,
    )
    bestObjective = bestBlendedScore
    cosineComponents = []
    selectedRowIndexes = []

    for _, row in paramDf.iloc[:maxFunctions].iterrows():
        candidateComponent = row["Amplitude"] * (1 + np.cos(np.radians(row["Period"] * angle - row["Phase"])))
        candidateSignal = reconstructedSignal + candidateComponent
        candidateTorsionScore, candidateTotalScore, candidateBlendedScore = _score_candidate(
            candidateSignal,
            qmTorsionEnergy,
            qmTotalEnergy=qmTotalEnergy,
            mmTotalEnergy=mmTotalEnergy,
            mmTorsionEnergy=mmTorsionEnergy,
            torsionObjectiveWeight=torsionObjectiveWeight,
            totalObjectiveWeight=totalObjectiveWeight,
        )
        scoreImprovement = bestBlendedScore - candidateBlendedScore
        if scoreImprovement < minScoreImprovement:
            continue

        candidateNTerms = len(cosineComponents) + 1
        candidateObjective = candidateBlendedScore + (complexityPenalty * (candidateNTerms / max(maxFunctions, 1)))
        if candidateObjective < bestObjective:
            reconstructedSignal = candidateSignal
            bestTorsionScore = candidateTorsionScore
            bestTotalScore = candidateTotalScore
            bestBlendedScore = candidateBlendedScore
            bestObjective = candidateObjective
            selectedRowIndexes.append(row.name)
            cosineComponents.append((len(cosineComponents) + 1, candidateComponent))
            if acceptableFitScore is not None:
                hasTotalScore = bestTotalScore is not None
                reachedTarget = (
                    bestTorsionScore <= acceptableFitScore and (not hasTotalScore or bestTotalScore <= acceptableFitScore)
                )
                if reachedTarget:
                    # stop at the first acceptable model to keep cosine count minimal
                    break

    selectedParamDf = paramDf.loc[selectedRowIndexes]
    return selectedParamDf, cosineComponents, len(cosineComponents)

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def calculate_rmsd(signal1: np.array, signal2: np.array) -> float:
    """Calculate the RMSD between two signals."""
    return np.sqrt(np.mean((signal1 - signal2) ** 2))
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def pad_energy_data(energyData: np.array, paddingFactor: int) -> np.array:
    """Tile energy data for Fourier padding."""
    return np.tile(energyData, paddingFactor)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def pad_angle_data(angleData: np.array, paddingFactor: int) -> np.array:
    """Tile angle data for Fourier padding."""
    return np.tile(angleData, paddingFactor)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rescale_torsion_angles(angle) -> float:
    """Wrap a torsion angle into the 0-360 degree range."""
    angle = angle % 360  # Reduce the angle to the 0-360 range
    if angle < 0:
        angle += 360  # Shift negative angles to the positive side
    return angle

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def perform_rfft(signal: np.array) -> np.array:
    """Perform a real FFT on the signal."""
    return np.fft.rfft(signal)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def get_frequencies(signalLength: int, sampleSpacing: int) -> np.array:
    """Return FFT frequency bins for the given signal length."""
    return np.fft.rfftfreq(signalLength, sampleSpacing)
##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def compute_amplitude_and_phase(fftResult: np.array, signalLength: int) -> Tuple[np.array, np.array]:
    """Convert FFT coefficients into amplitude and phase arrays."""
    amplitudes: np.array = np.abs(fftResult) * 2 / signalLength
    phases: np.array = np.angle(fftResult)

    return amplitudes, phases

##🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    
    raise NotImplementedError

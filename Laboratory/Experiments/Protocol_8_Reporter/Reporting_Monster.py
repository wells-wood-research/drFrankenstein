import pandas as pd
import os.path as p
from copy import deepcopy
from shutil import copy
import os

from . import Reporting_Assistant

def process_fitting_results(config: dict) -> dict:
    ## unpack config
    imagesDir = config["runtimeInfo"]["madeByReporting"]["imagesDir"]
    qmmmParameterFittingDir = config["runtimeInfo"]["madeByStitching"]["qmmmParameterFittingDir"]
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    # nShuffles = config["parameterFittingInfo"]["nShuffles"]
    shufflesCompleted = config["runtimeInfo"]["madeByStitching"]["shufflesCompleted"]
    ## make directories
    fittingImagesDir = p.join(imagesDir, "fitting_images")    
    os.makedirs(fittingImagesDir, exist_ok=True)


    fittingData = {
        "nShuffles": shufflesCompleted,
        "torsionsToScan": config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"],
        "finalParameters": config["runtimeInfo"]["madeByStitching"]["finalParameters"],
        "maxCosineFunctions": config["parameterFittingInfo"]["maxCosineFunctions"],
        "sagvolSmoothing": config["parameterFittingInfo"]["sagvolSmoothing"],
        "l2DampingFactor": config["parameterFittingInfo"]["l2DampingFactor"],
    }

    ## collect all torsions mae png
    allTorsionMaePng = p.join(qmmmParameterFittingDir, "run_mean_average_error.png")
    destAllMaePng = p.join(fittingImagesDir, "all_torsions_mae.png")
    copy(allTorsionMaePng, destAllMaePng)
    relativeAllMaePng = p.relpath(destAllMaePng, reporterDir)
    fittingData["allTorsionMaePng"] = relativeAllMaePng

    ## collect pngs and gifs
    fittingImages = {}
    for torsionTag in os.listdir(qmmmParameterFittingDir):
        torsionFittingDir = p.join(qmmmParameterFittingDir, torsionTag)
        if not p.isdir(torsionFittingDir):
            continue
        maePng = p.join(torsionFittingDir, f"mean_average_error.png")
        finalPng = p.join(torsionFittingDir, f"fitting_shuffle_{shufflesCompleted}.png")
        fittingGif = p.join(torsionFittingDir, f"torsion_fitting.gif")
        destPng = p.join(fittingImagesDir, f"{torsionTag}_final.png")
        destGif = p.join(fittingImagesDir, f"{torsionTag}_fitting.gif")
        destMaePng = p.join(fittingImagesDir, f"{torsionTag}_mae.png")
        copy(finalPng, destPng)
        copy(fittingGif, destGif)
        copy(maePng, destMaePng)
        relativePng = p.relpath(destPng, reporterDir)
        relativeGif = p.relpath(destGif, reporterDir)
        relativeMaePng = p.relpath(destMaePng, reporterDir)
        fittingImages[torsionTag] = {
            "finalPng": relativePng,
            "fittingGif": relativeGif,
            "maePng": relativeMaePng,
        }
    fittingData["fittingImages"] = fittingImages
    return fittingData




def process_charges_results(config: dict) -> dict:
    ## unpack config
    imagesDir = config["runtimeInfo"]["madeByReporting"]["imagesDir"]
    chargesImagesDir = p.join(imagesDir, "charges_images")    
    os.makedirs(chargesImagesDir, exist_ok=True)

    chargesData = {
        "chargeFittingProtocol": config["chargeFittingInfo"]["chargeFittingProtocol"],
        "optMethod" : config["chargeFittingInfo"]["optMethod"],
        "optSolvationMethod" : config["chargeFittingInfo"]["optSolvationMethod"],   
        "singlePointMethod" : config["chargeFittingInfo"]["singlePointMethod"], 
        "singlePointSolvationMethod" : config["chargeFittingInfo"]["singlePointSolvationMethod"],
    }
     
    if config["chargeFittingInfo"]["nConformers"] == -1:
        chargesData["nConformers"] = config["runtimeInfo"]["madeByConformers"]["nConformersGenerated"]
    else:
        chargesData["nConformers"] = config["chargeFittingInfo"]["nConformers"]

    chargesCsv = config["runtimeInfo"]["madeByCharges"]["chargesCsv"]
    chargesDf = pd.read_csv(chargesCsv)

    partialChargeData = chargesDf.loc[:, ["ATOM_NAME", "Charge"]].values.tolist()
    partialChargeHeaders = ["Atom Name", "Partial Charge"]

    chargesData["partialChargeData"] = partialChargeData
    chargesData["partialChargeHeaders"] = partialChargeHeaders

    chargeHtml, minCharge, maxCharge = Reporting_Assistant.make_charge_visualisation(config, chargesImagesDir)

    colorBar = Reporting_Assistant.make_charge_color_bar(minCharge, maxCharge, chargesImagesDir, config)

    chargesData["colorBar"] = colorBar
    chargesData["chargeHtml"] = chargeHtml

    return chargesData





def process_twist_results(config: dict) -> dict:
    # Unpack config
    imagesDir = config["runtimeInfo"]["madeByReporting"]["imagesDir"]
    torsionImagesDir = p.join(imagesDir, "torsion_images")
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    os.makedirs(torsionImagesDir, exist_ok=True)

    twistData = {
        "scanMethod": config["torsionScanInfo"]["scanMethod"],
        "scanSolvationMethod": config["torsionScanInfo"]["scanSolvationMethod"],
        "singlePointMethod": config["torsionScanInfo"]["singlePointMethod"],
        "singlePointSolvationMethod": config["torsionScanInfo"]["singlePointSolvationMethod"],
        "nRotatableBonds": config["runtimeInfo"]["madeByTwisting"]["nRotatableBonds"],
        "torsionsToScan": deepcopy(config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"])  # Deep copy
    }
    if config["torsionScanInfo"]["nConformers"] == -1:
        twistData["nConformers"] = config["runtimeInfo"]["madeByConformers"]["nConformersGenerated"]
    else:
        twistData["nConformers"] = config["torsionScanInfo"]["nConformers"]

    torsionsToScan = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]

    for torsionTag, torsionData in torsionsToScan.items():
        for pngTag in ["scanPng", "spPng", "scanVsSpPng"]:
            if torsionData[pngTag] is None:
                continue
            else:
                pngFile = torsionData[pngTag]
                destPng = p.join(torsionImagesDir, f"{torsionTag}_{pngTag}.png")
                copy(pngFile, destPng)  
                relativePath = p.relpath(destPng, reporterDir)
                twistData["torsionsToScan"][torsionTag][pngTag] = relativePath

    torsionHtmls = Reporting_Assistant.make_highlighted_torsion_visualisations(config, torsionImagesDir)
    twistData["torsionHtmls"] = torsionHtmls    

    return twistData



def process_wriggle_results(config: dict) -> dict:
    ## unpack config
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    imagesDir = config["runtimeInfo"]["madeByReporting"]["imagesDir"]
    conformerImagesDir = p.join(imagesDir, "conformer_images")
    os.makedirs(conformerImagesDir, exist_ok=True)


    conformerHtmls = Reporting_Assistant.make_conformer_visualisations(config, conformerImagesDir)

    wriggleData = {
        "conformerHtmls": conformerHtmls,
        "conformerEnergies": config["runtimeInfo"]["madeByConformers"]["conformerEnergies"],
        "nConformersGenerated": config["runtimeInfo"]["madeByConformers"]["nConformersGenerated"]
    }

    return wriggleData

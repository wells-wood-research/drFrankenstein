import pandas as pd
import py3Dmol
from pdbUtils import pdbUtils # Your custom PDB parsing utility
import os
from os import path as p
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

def make_charge_visualisation(config, outDir):
    ## unpack config
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    chargesCsv = config["runtimeInfo"]["madeByCharges"]["chargesCsv"]
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    cappedDf = pdbUtils.pdb2df(cappedPdb)
    chargesDf = pd.read_csv(chargesCsv, index_col="Unnamed: 0")

    chargeGradientColors = make_charge_gradient_colours()

    cappedDf["BETAFACTOR"] = chargesDf["Charge"]

    annotatedPdb = p.join(outDir, "annotated_charged.pdb")
    pdbUtils.df2pdb(cappedDf, annotatedPdb)

    minCharge = chargesDf["Charge"].min()
    maxCharge = chargesDf["Charge"].max()

    ## load PDB into a string
    with open(annotatedPdb, 'r') as f:
        pdbBlock = "".join(f.readlines())


    # Create 3D visualization with py3Dmol
    view = py3Dmol.view(width=600, height=600) # Increased size slightly for better viewing
    view.addModel(pdbBlock, 'pdb', {"keepH": True})
    view.setStyle({
        'stick': {
            'radius': 0.1,
            'color': 'green'
        },
        'sphere': {
            'scale': 0.5,
            "opacity": 0.75,
            'colorscheme': {'prop': 'b', 'gradient': "linear", "colors": chargeGradientColors }
        }
    })
    view.setBackgroundColor("black")

    for index, row in chargesDf.iterrows():
        view.addLabel(row["ATOM_NAME"],
                        {'fontSize': 10, 'fontColor': 'yellow', 
                        'showBackground': True, 'backgroundOpacity': 1, 'backgroundColor': 'black',
                        'alignment': 'topCenter'},
                        {'serial': index+1}) 
        view.addLabel(f"{round(row['Charge'],2):.2f}",
                        {'fontSize': 10, 'fontColor': 'yellow', 
                        'showBackground': True, 'backgroundOpacity': 1, 'backgroundColor': 'black',
                        'alignment': 'bottomCenter'},
                        {'serial': index+1}) 
    view.zoomTo()

    # Save the visualization as HTML
    chargeHtml = p.join(outDir, "charges.html")
    with open(chargeHtml, 'w') as f:
        f.write(view._make_html())

    os.remove(annotatedPdb)

    relativePath = p.relpath(chargeHtml, reporterDir)

    return relativePath, minCharge, maxCharge


def make_charge_color_bar(minCharge, maxCharge, outDir, config):
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    chargeGradientColors = make_charge_gradient_colours() # Assuming this function is defined

    # Create figure with black background
    fig, ax = plt.subplots(figsize=(6, 1), facecolor='black') # Figure size in inches
    fig.subplots_adjust(bottom=0.5) # Adjust subplot to leave space for label below color bar

    cmap = mcolors.ListedColormap(chargeGradientColors)
    norm = mcolors.Normalize(vmin=minCharge, vmax=maxCharge)

    # Create ScalarMappable
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    # Create colorbar with white elements
    cb = fig.colorbar(mappable, cax=ax, orientation='horizontal')
    
    # Set label and customize appearance
    cb.set_label('Charge', color='white', fontsize=12)
    
    # Customize colorbar appearance
    cb.outline.set_edgecolor('white')  # White outline
    cb.ax.tick_params(colors='white')  # White tick marks and labels
    ax.set_facecolor('black')  # Ensure axes background is black

    # Set figure background to black (redundant but ensures consistency)
    fig.set_facecolor('black')

    outputFilePath = os.path.join(outDir, "charge_color_bar.png")
    
    # Save the figure with a specific DPI to control pixel height
    plt.savefig(outputFilePath, facecolor='black', edgecolor='none')
    plt.close(fig) # Close the figure to free memory

    relativePath = p.relpath(outputFilePath, reporterDir)
    return relativePath


def  make_conformer_visualisations(config, outDir):
    ## unpack config
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    reporterDir =     config["runtimeInfo"]["madeByReporting"]["reporterDir"]

    vibrantColors = make_vibrant_colors()

    conformerHtmls = []
    for idx, xyzFile in enumerate(conformerXyzs):
        name = xyzFile.split("/")[-1].split(".")[0]
        xyzBlock = xyz2String(xyzFile)
        # Create 3D visualization with py3Dmol
        view = py3Dmol.view(width=200, height=200) 
        view.addModel(xyzBlock, 'xyz', {"keepH": True})
        # Use modulo to cycle through colors
        color = vibrantColors[idx % len(vibrantColors)]
        view.setStyle({"elem": "C"},{
            'stick': {
                'radius': 0.2, "color": color
            },
            'sphere': {
                'scale': 0.3, "color": color
            }
        })
        view.setStyle({"elem": "C", "invert": True},{
            'stick': {
                'radius': 0.2
            },
            'sphere': {
                'scale': 0.3
            }
        })

        view.setBackgroundColor("black")

        # Save the visualization as HTML
        conformerHtml = p.join(outDir, f"{name}.html")
        with open(conformerHtml, 'w') as f:
            f.write(view._make_html())

        relativePath = p.relpath(conformerHtml, reporterDir)

        conformerHtmls.append(relativePath)

    return conformerHtmls


def make_highlighted_torsion_visualisations(config, outDir):
    ## unpack config
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    uniqueRotatableDihedrals = config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"]
    reporterDir =     config["runtimeInfo"]["madeByReporting"]["reporterDir"]

    pdbBlock = pdb2string(cappedPdb)

    htmlFiles = {}
    for torsionTag, dihedralData in uniqueRotatableDihedrals.items():
        view = py3Dmol.view(width=300, height=300)
        view.addModel(pdbBlock, "pdb",  {"keepH": True})

        ## simple whiteCarbons balls and sticks
        view.setStyle({},{
            'stick': {
                'radius': 0.2
            },
            'sphere': {
                'scale': 0.3
            }
        })
        view.setStyle({"atom": dihedralData["ATOM_NAMES"]}, {
            'sphere': {
                'scale': 0.3,
                'color': "#00FF00"
            },
                    'stick': {
                'radius': 0.2
            }
        })
        view.setBackgroundColor("black")

        # Save the visualization as HTML
        htmlFile = p.join(outDir, f"{torsionTag}.html")
        with open(htmlFile, 'w') as f:
            f.write(view._make_html())

        relativePath = p.relpath(htmlFile, reporterDir)
        htmlFiles[torsionTag] = relativePath

    return htmlFiles

def make_vibrant_colors():
    vibrantColours = [
    "#FF0000",  # Bright Red
    "#00FF00",  # Lime Green
    "#0000FF",  # Pure Blue
    "#FFFF00",  # Bright Yellow
    "#FF00FF",  # Magenta
    "#00FFFF",  # Cyan
    "#FF4500",  # Orange Red
    "#32CD32",  # Lime Green
    "#1E90FF",  # Dodger Blue
    "#FFD700",  # Gold
    "#FF69B4",  # Hot Pink
    "#00CED1",  # Dark Turquoise
    "#FF6347",  # Tomato
    "#ADFF2F",  # Green Yellow
    "#4682B4",  # Steel Blue
    "#FFA500",  # Orange
    "#DA70D6",  # Orchid
    "#20B2AA",  # Light Sea Green
    "#DC143C",  # Crimson
    "#7CFC00"   # Lawn Green
    ]
    return vibrantColours


def make_charge_gradient_colours():
    chargeGradientColors = ['#FF00FF', '#FF66FF', '#FFCCFF', '#FFF5FF', '#FFFFFF', '#F5FFF5', '#CCFFCC', '#33FF33', '#00FF00']
    return chargeGradientColors

def pdb2string(pdbFile):
    with open(pdbFile, "r") as f:
        pdbBlock = "".join(f.readlines())
    return pdbBlock


def xyz2String(xyzFile):
    with open(xyzFile, 'r') as f:
        xyzBlock = "".join(f.readlines())
    return xyzBlock

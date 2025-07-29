import pandas as pd
import py3Dmol
from pdbUtils import pdbUtils # Your custom PDB parsing utility
from shutil import copy
import os
from os import path as p
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import matplotlib.text as mtext # For creating text artists
from matplotlib.legend_handler import HandlerBase # Base class for custom handlers
import matplotlib.ticker as ticker


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as mpe
import numpy as np
import seaborn as sns
import os
import os.path as p

# Utility Functions
def format_time_hms(seconds):
    """Converts seconds to HH:MM:SS string format."""
    seconds = int(round(seconds))
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    secs = seconds % 60
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"

def pdb_to_string(pdbFile):
    with open(pdbFile, "r") as f:
        pdbBlock = "".join(f.readlines())
    return pdbBlock

def xyz_to_string(xyzFile):
    with open(xyzFile, 'r') as f:
        xyzBlock = "".join(f.readlines())
    return xyzBlock

def make_vibrant_colors():
    vibrantColors = [
        "#FF0000",  # Bright Red
        "#00FF00",  # Lime Green
        "#0000FF",  # Pure Blue
        "#FFFF00",  # Bright Yellow
        "#FF00FF",  # Magenta
        "#00FFFF",  # Cyan
        "#FF4500",  # Orange Red
        "#32CD32",  # Lime Green (repeated, consider variety if an issue)
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
    return vibrantColors

def make_charge_gradient_colors():
    chargeGradientColors = ['#FF00FF', '#FF66FF', '#FFCCFF', '#FFF5FF', '#FFFFFF', '#F5FFF5', '#CCFFCC', '#33FF33', '#00FF00']
    return chargeGradientColors

# Gantt Chart Helper Functions
def initialize_matplotlib_settings(defaultFontSize):
    """Sets global Matplotlib rcParams and verifies font."""
    plt.rcParams['font.family'] = 'monospace'
    plt.rcParams['font.monospace'] = ['Consolas', 'DejaVu Sans Mono', 'Courier New']
    plt.rcParams['font.size'] = defaultFontSize

def create_figure_and_axes(numFunctions, backgroundColor, minFigHeight, figHeightPerFunction, figHeightBasePadding):
    """Creates the Matplotlib figure and axes with basic styling."""
    figHeight = max(minFigHeight, numFunctions * figHeightPerFunction + figHeightBasePadding)
    fig, ax = plt.subplots(figsize=(24, figHeight))
    fig.patch.set_facecolor(backgroundColor)
    ax.set_facecolor(backgroundColor)
    return fig, ax

def draw_task_bars(ax, timeDf, barColors, totalRuntimeSeconds, textColor, backgroundColor, elementColor, defaultFontSize):
    """Draws the Gantt chart bars for each task."""
    numFunctions = len(timeDf)
    for i, row in timeDf.iterrows():
        yPos = numFunctions - 1 - i
        barStartMinutes = row['startTimeSeconds'] / 60
        barDurationMinutes = row['executionTimeSeconds'] / 60
        
        ax.broken_barh(
            xranges=[(barStartMinutes, barDurationMinutes)],
            yrange=(yPos - 0.4, 0.8),
            facecolors=[barColors[i % len(barColors)] if numFunctions > 0 else 'blue'],
            edgecolor=elementColor,
            linewidth=1
        )
        if totalRuntimeSeconds > 0 and row['executionTimeSeconds'] > (0.05 * totalRuntimeSeconds):
            timeLabel = format_time_hms(row['executionTimeSeconds'])
            textXPosMinutes = barStartMinutes + barDurationMinutes / 2
            ax.text(textXPosMinutes,
                    yPos, timeLabel,
                    ha='center', va='center', color=textColor, fontsize=defaultFontSize,
                    path_effects=[mpe.Stroke(linewidth=2, foreground=backgroundColor), mpe.Normal()])

def draw_group_annotations(ax, timeDf, functionMap, numFunctions,
                           groupLabelColor, groupBoxEdgeColor,
                           groupBoxPaddingY, groupLabelYOffset, 
                           groupLabelMaxCharsPerLine, 
                           estimatedLineHeightForGroupLabel, 
                           groupLabelFontSize):
    """Draws group boxes and numerical labels (in a styled box) on the Gantt chart."""
    yMaxFromGroupLabels = 0.0
    groupLegendDetails = [] 

    if numFunctions == 0:
        yMaxFromGroupLabels = 1.0

    for groupIdx, (groupName, functionsInGroupMap) in enumerate(functionMap.items()):
        csvNamesInGroup = list(functionsInGroupMap.keys())
        groupMemberRows = timeDf[timeDf['functionName'].isin(csvNamesInGroup)]
        if groupMemberRows.empty:
            continue

        memberIndices = groupMemberRows.index.tolist()
        if not memberIndices:
            continue

        currentMinIndexInDf = min(memberIndices)
        currentMaxIndexInDf = max(memberIndices)

        yPosOfTopBarInGroup = numFunctions - 1 - currentMinIndexInDf
        yPosOfBottomBarInGroup = numFunctions - 1 - currentMaxIndexInDf

        boxBottomYCoord = yPosOfBottomBarInGroup - 0.4 - groupBoxPaddingY
        boxTopYCoord = yPosOfTopBarInGroup + 0.4 + groupBoxPaddingY

        rectYStart = boxBottomYCoord
        rectHeight = boxTopYCoord - rectYStart
        if rectHeight <= 0:
            continue

        boxXStartSeconds = groupMemberRows['startTimeSeconds'].min()
        boxXEndSeconds = groupMemberRows['endTimeSeconds'].max()

        rectXStartMinutes = boxXStartSeconds / 60
        rectWidthMinutes = (boxXEndSeconds - boxXStartSeconds) / 60
        if rectWidthMinutes < 0:
            rectWidthMinutes = 0

        rect = patches.Rectangle(
            (rectXStartMinutes, rectYStart),
            rectWidthMinutes,
            rectHeight,
            linewidth=0.75,
            edgecolor=groupBoxEdgeColor,
            facecolor='none',
            zorder=0.5
        )
        ax.add_patch(rect)

        labelXPos = rectXStartMinutes + rectWidthMinutes + 0.1
        labelYPos = rectYStart + rectHeight / 2

        groupIdToDisplay = str(groupIdx + 1) 
        finalLabelTextForGroup = groupIdToDisplay

        groupLegendDetails.append((groupIdToDisplay, groupName)) 

        bboxProps = dict(boxstyle="square,pad=0.2", 
                          facecolor="black",          
                          edgecolor="lime",           
                          linewidth=1)                

        ax.text(labelXPos,
                labelYPos,
                finalLabelTextForGroup, 
                ha='left',              
                va='center',            
                color=groupLabelColor,  
                fontsize=groupLabelFontSize,
                weight='bold',
                clip_on=False,          
                bbox=bboxProps         
                )
        yMaxFromGroupLabels = max(yMaxFromGroupLabels, boxTopYCoord + 0.1)

    return yMaxFromGroupLabels, groupLegendDetails

class HandlerStyledIdBox(HandlerBase):
    """Custom legend handler to draw a text ID within a styled box."""
    def __init__(self, textColorName, boxFacecolorName, boxEdgecolorName,
                fontWeight='bold', fontFamily='sans-serif', 
                boxStyleStr='square,pad=0.3', 
                fontSizeFactor=0.9, 
                boxLinewidth=1):
        super().__init__()
        self.textColorName = textColorName
        self.boxFacecolorName = boxFacecolorName
        self.boxEdgecolorName = boxEdgecolorName
        self.fontWeight = fontWeight
        self.fontFamily = fontFamily
        self.boxStyleStr = boxStyleStr
        self.fontSizeFactor = fontSizeFactor
        self.boxLinewidth = boxLinewidth

    def create_artists(self, legend, origHandle, 
                    xdescent, ydescent, width, height, fontsize, trans):
        
        groupIdStr = str(origHandle) 

        textArtist = mtext.Text(
            width / 2.0, height / 2.0, 
            groupIdStr,
            color=self.textColorName,
            fontsize=fontsize * self.fontSizeFactor, 
            ha="center",                             
            va="center",                             
            fontweight=self.fontWeight,
            fontfamily=self.fontFamily,
            transform=trans 
        )

        bboxProps = dict(
            boxstyle=self.boxStyleStr,
            facecolor=self.boxFacecolorName,
            edgecolor=self.boxEdgecolorName,
            linewidth=self.boxLinewidth
        )
        textArtist.set_bbox(bboxProps)
        
        return [textArtist]

def finalize_plot_styling(fig, ax, timeDf, displayNameMap, yMaxForPlot, 
                          textColor, elementColor, backgroundColor, defaultFontSize,
                          groupLegendDetails, groupLabelColor, groupBoxEdgeColor):
    """Applies final styling to the plot (labels, ticks, limits, grid, spines, and group legend)."""
    numFunctions = len(timeDf)

    ax.set_title('Function Execution Gantt Chart', fontsize=24, pad=20, weight='bold', color=textColor)

    if numFunctions > 0:
        ax.set_yticks(np.arange(numFunctions))
        originalYTickLabels = timeDf['functionName'].iloc[::-1].tolist() 
        displayYTickLabels = [displayNameMap.get(name, name) for name in originalYTickLabels]
        ax.set_yticklabels(displayYTickLabels, fontsize=defaultFontSize)
    else:
        ax.set_yticks([])

    maxEndTimeValSeconds = 0.0
    if not timeDf.empty:
        maxEndTimeValSeconds = timeDf['endTimeSeconds'].max()
        if pd.isna(maxEndTimeValSeconds): 
            maxEndTimeValSeconds = 0.0
    ## decide if xTicks will be in hours or minutes
    xLimUpperMinutes = 1.0 

    time_data_seconds = np.linspace(0, maxEndTimeValSeconds, 200)
    if maxEndTimeValSeconds < 18000:  # 5 hours (5 * 3600 = 18000 seconds)
        print("Plotting: less than 5 hours, using minutes")
        # Convert data and limits to minutes
        xLimUpperDisplay = (maxEndTimeValSeconds / 60.0) * 1.05

        ax.set_xlim(0, xLimUpperDisplay)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=30))  # Ticks every 30 minutes
        ax.set_xlabel('Time (minutes)', fontsize=defaultFontSize, labelpad=15, color=textColor)

    else:  # >= 5 hours
        print("Plotting: more than 5 hours, using hours")
        # Convert data and limits to hours
        xLimUpperDisplay = (maxEndTimeValSeconds / 3600.0) * 1.05

        ax.set_xlim(0, xLimUpperDisplay)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1))  # Ticks every 1 hour
        ax.set_xlabel('Time (hours)', fontsize=defaultFontSize, labelpad=15, color=textColor)


    currentYLimBottom, _ = ax.get_ylim() 
    ax.set_ylim(currentYLimBottom, yMaxForPlot)

    ax.grid(True, axis='x', linestyle=':', alpha=0.6, color=elementColor)
    ax.grid(False, axis='y') 
    ax.tick_params(axis='x', colors=textColor, labelsize=defaultFontSize)
    ax.tick_params(axis='y', colors=textColor, labelsize=defaultFontSize) 
    ax.spines['bottom'].set_color(elementColor)
    ax.spines['left'].set_color(elementColor)
    ax.spines['top'].set_color(backgroundColor) 
    ax.spines['right'].set_color(backgroundColor) 
    sns.despine(top=True, right=True, left=False, bottom=False)

    if groupLegendDetails: 
        legendHandles = []
        legendLabels = []
        
        styledIdKeyHandler = HandlerStyledIdBox(
            textColorName=groupLabelColor,    
            boxFacecolorName='black',         
            boxEdgecolorName=groupLabelColor, 
            fontSizeFactor=0.85,              
            boxStyleStr='square,pad=0.3',     
            boxLinewidth=1
        )
        
        handlerMapForLegend = {}

        for groupIdStr, fullGroupName in groupLegendDetails:
            legendHandles.append(groupIdStr)
            legendLabels.append(f"{fullGroupName}") 
            handlerMapForLegend[groupIdStr] = styledIdKeyHandler

        ax.legend(
            legendHandles,  
            legendLabels,   
            fontsize=defaultFontSize, 
            title_fontsize=defaultFontSize * 0.9,
            loc='lower left',              
            bbox_to_anchor=(0.0, 0.0),     
            facecolor=backgroundColor,     
            edgecolor="lime",        
            labelcolor=textColor,          
            handler_map=handlerMapForLegend 
        )
    fig.subplots_adjust(left=0.22, right=0.95, top=0.95, bottom=0.15)

# Main Visualization Orchestration Functions
def generate_gantt_chart(config):
    """
    Generates and saves the Gantt chart based on provided function data dictionary
    and saves it to a path relative to reporterDir.
    """
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    timeInfo  = config["runtimeInfo"]["timeInfo"]
    
    backgroundColor = 'black'
    textColor = 'yellow'
    elementColor = 'yellow'
    groupLabelColor = 'lime'
    groupBoxEdgeColor = 'lime'
    defaultFontSize = 16
    groupLabelFontSize = 16
    groupLabelMaxCharsPerLine = 25 
    estimatedLineHeightForGroupLabel = 0.4
    groupBoxPaddingY = 0.10
    groupLabelYOffset = 0.1
    minFigHeight = 9
    figHeightPerFunction = 0.5
    figHeightBasePadding = 2

    dataForDf = []
    for funcNameKey, detailsDict in timeInfo.items():
        dataForDf.append({
            'functionName': funcNameKey,
            'executionTimeSeconds': detailsDict['executionTime']
        })
    timeDf = pd.DataFrame(dataForDf)

    if timeDf.empty:
        timeDf = pd.DataFrame(columns=['functionName', 'executionTimeSeconds', 'startTimeSeconds', 'endTimeSeconds'])
        totalRuntimeSeconds = 0.0
    else:
        timeDf['executionTimeSeconds'] = pd.to_numeric(timeDf['executionTimeSeconds'], errors='coerce').fillna(0)
        timeDf['startTimeSeconds'] = timeDf['executionTimeSeconds'].cumsum().shift(1).fillna(0)
        timeDf['endTimeSeconds'] = timeDf['startTimeSeconds'] + timeDf['executionTimeSeconds']
        totalRuntimeSeconds = timeDf['executionTimeSeconds'].sum()

    displayNameMap = {}
    functionMap = {} 
    for funcNameKey, detailsDict in timeInfo.items():
        csvName = funcNameKey
        displayName = detailsDict['functionAlias']
        groupName = detailsDict['functionGroup']
        displayNameMap[csvName] = displayName
        if groupName not in functionMap:
            functionMap[groupName] = {}
        functionMap[groupName][csvName] = displayName

    timeImagesDir = p.join(reporterDir, "Images", "time_data")
    os.makedirs(timeImagesDir, exist_ok=True)

    initialize_matplotlib_settings(defaultFontSize)
    numFunctions = len(timeDf)
    fig, ax = create_figure_and_axes(numFunctions, backgroundColor, minFigHeight, figHeightPerFunction, figHeightBasePadding)
    
    barColors = sns.color_palette('viridis', n_colors=max(1, numFunctions))

    draw_task_bars(ax, timeDf, barColors, totalRuntimeSeconds, textColor, backgroundColor, elementColor, defaultFontSize)
    
    yMaxForPlot = numFunctions if numFunctions > 0 else 1.0
    yMaxFromGroupLabels, groupLegendDetails = draw_group_annotations(ax, timeDf, functionMap, numFunctions, 
                                               groupLabelColor, groupBoxEdgeColor, 
                                               groupBoxPaddingY, groupLabelYOffset, 
                                               groupLabelMaxCharsPerLine, 
                                               estimatedLineHeightForGroupLabel, 
                                               groupLabelFontSize)
    yMaxForPlot = max(yMaxForPlot, yMaxFromGroupLabels)
    
    finalize_plot_styling(fig, ax, timeDf, displayNameMap, yMaxForPlot, 
                          textColor, elementColor, backgroundColor, defaultFontSize,
                          groupLegendDetails, groupLabelColor, groupBoxEdgeColor)

    ganttPngPath = p.join(timeImagesDir, "gantt_chart.png")
    plt.savefig(ganttPngPath, dpi=300, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close(fig)

    relativePath = p.relpath(ganttPngPath, reporterDir)
    return relativePath

def make_charge_visualisation(config, outDir):
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    chargesCsv = config["runtimeInfo"]["madeByCharges"]["chargesCsv"]
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    
    cappedDf = pdbUtils.pdb2df(cappedPdb)
    chargesDf = pd.read_csv(chargesCsv, index_col="Unnamed: 0")

    chargeGradientColors = make_charge_gradient_colors()

    cappedDf["BETAFACTOR"] = chargesDf["Charge"]

    annotatedPdb = p.join(outDir, "annotated_charged.pdb")
    pdbUtils.df2pdb(cappedDf, annotatedPdb)

    minCharge = chargesDf["Charge"].min()
    maxCharge = chargesDf["Charge"].max()

    pdbBlock = pdb_to_string(annotatedPdb)

    view = py3Dmol.view(width=600, height=600) 
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

    chargeHtml = p.join(outDir, "charges.html")
    with open(chargeHtml, 'w') as f:
        f.write(view._make_html())

    os.remove(annotatedPdb)
    relativePath = p.relpath(chargeHtml, reporterDir)
    return relativePath, minCharge, maxCharge

def make_charge_color_bar(minCharge, maxCharge, outDir, config):
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    chargeGradientColors = make_charge_gradient_colors() 

    fig, ax = plt.subplots(figsize=(6, 1), facecolor='black') 
    fig.subplots_adjust(bottom=0.5) 

    cmap = mcolors.ListedColormap(chargeGradientColors)
    norm = mcolors.Normalize(vmin=minCharge, vmax=maxCharge)

    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    cb = fig.colorbar(mappable, cax=ax, orientation='horizontal')
    
    cb.set_label('Charge', color='white', fontsize=12)
    
    cb.outline.set_edgecolor('white')  
    cb.ax.tick_params(colors='white')  
    ax.set_facecolor('black')  

    fig.set_facecolor('black')

    outputFilePath = os.path.join(outDir, "charge_color_bar.png")
    
    plt.savefig(outputFilePath, facecolor='black', edgecolor='none', bbox_inches='tight')
    plt.close(fig) 

    relativePath = p.relpath(outputFilePath, reporterDir)
    return relativePath

def make_conformer_visualisations(config, outDir):
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    reporterDir =     config["runtimeInfo"]["madeByReporting"]["reporterDir"]

    vibrantColors = make_vibrant_colors()

    conformerHtmls = []
    for idx, xyzFile in enumerate(conformerXyzs):
        name = xyzFile.split("/")[-1].split(".")[0]
        xyzBlock = xyz_to_string(xyzFile)
        view = py3Dmol.view(width=200, height=200) 
        view.addModel(xyzBlock, 'xyz', {"keepH": True})
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

        conformerHtml = p.join(outDir, f"{name}.html")
        with open(conformerHtml, 'w') as f:
            f.write(view._make_html())

        relativePath = p.relpath(conformerHtml, reporterDir)
        conformerHtmls.append(relativePath)

    return conformerHtmls

def make_highlighted_torsion_visualisations(config, outDir):
    cappedPdb = config["runtimeInfo"]["madeByCapping"]["cappedPdb"]
    uniqueRotatableDihedrals = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]
    reporterDir =     config["runtimeInfo"]["madeByReporting"]["reporterDir"]

    pdbBlock = pdb_to_string(cappedPdb)

    htmlFiles = {}
    for torsionTag, dihedralData in uniqueRotatableDihedrals.items():
        view = py3Dmol.view(width=600, height=300)
        view.addModel(pdbBlock, "pdb",  {"keepH": True})

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
                'color': "#00FF00" # Lime Green
            },
                    'stick': {
                'radius': 0.2
            }
        })

        for atomName, atomIndex in zip(dihedralData["ATOM_NAMES"], dihedralData["ATOM_INDEXES"]):
            view.addLabel(atomName,
                            {'fontSize': 10, 'fontColor': 'yellow', 
                            'showBackground': True, 'backgroundOpacity': 1, 'backgroundColor': 'black',
                            'alignment': 'center'},
                            {'serial': atomIndex + 1}) 

        view.setBackgroundColor("black")

        htmlFile = p.join(outDir, f"{torsionTag}.html")
        with open(htmlFile, 'w') as f:
            f.write(view._make_html())

        relativePath = p.relpath(htmlFile, reporterDir)
        htmlFiles[torsionTag] = relativePath

    return htmlFiles

def copy_images(config):
    imagesDir = config["runtimeInfo"]["madeByReporting"]["imagesDir"]
    
    thisDir = os.path.dirname(os.path.abspath(__file__))
    
    lightningJpg = p.join(thisDir, "templates", "Lightning.jpg")
    lightningJpgDest = p.join(imagesDir, "Lightning.jpg")
    copy(lightningJpg, lightningJpgDest)

    cobbleJpg = p.join(thisDir, "templates", "Gothic_Cobble.jpg")
    cobbleJpgDest = p.join(imagesDir, "Gothic_Cobble.jpg")
    copy(cobbleJpg, cobbleJpgDest)
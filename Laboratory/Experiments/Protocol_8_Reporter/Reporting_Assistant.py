import pandas as pd
import py3Dmol
from pdbUtils import pdbUtils # Your custom PDB parsing utility
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
import textwrap
from matplotlib.font_manager import findfont, FontProperties

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as mpe
import numpy as np
import seaborn as sns
import textwrap
import os
import os.path as p

def format_time_hms(seconds):
    """Converts seconds to HH:MM:SS string format."""
    seconds = int(round(seconds))
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    secs = seconds % 60
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"

def initialize_matplotlib_settings(defaultFontSize):
    """Sets global Matplotlib rcParams and verifies font."""
    plt.rcParams['font.family'] = 'monospace'
    plt.rcParams['font.monospace'] = ['Consolas', 'DejaVu Sans Mono', 'Courier New']
    plt.rcParams['font.size'] = defaultFontSize

# load_and_prepare_data is removed as its functionality for the new input
# format is integrated into generate_gantt_chart.

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
            linewidth=0.75
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
                           groupLabelMaxCharsPerLine, # Parameter kept for signature
                           estimatedLineHeightForGroupLabel, 
                           groupLabelFontSize):
    """Draws group boxes and numerical labels (in a styled box) on the Gantt chart."""
    yMaxFromGroupLabels = 0.0
    group_legend_details = [] # To store (group_id_str, full_group_name)

    if numFunctions == 0: 
        yMaxFromGroupLabels = 1.0 

    for group_idx, (groupName, functionsInGroupMap) in enumerate(functionMap.items()):
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
        if rectHeight <= 0: continue

        boxXStartSeconds = groupMemberRows['startTimeSeconds'].min()
        boxXEndSeconds = groupMemberRows['endTimeSeconds'].max()
        
        rectXStartMinutes = boxXStartSeconds / 60
        rectWidthMinutes = (boxXEndSeconds - boxXStartSeconds) / 60
        if rectWidthMinutes < 0: rectWidthMinutes = 0 

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

        labelXCenterMinutes = rectXStartMinutes + rectWidthMinutes / 2
        labelYPosBottomOfText = boxTopYCoord + groupLabelYOffset
        
        group_id_to_display = str(group_idx + 1) # Use numerical ID
        finalLabelTextForGroup = group_id_to_display
        
        group_legend_details.append((group_id_to_display, groupName)) # Store for legend
        
        numLinesInLabel = 1 # Numerical ID is single line
        
        # Define the bounding box properties for the numerical ID
        bbox_props = dict(boxstyle="square,pad=0.2", # Square box with some padding
                          facecolor="black",          # Black background
                          edgecolor="lime",           # Lime border
                          linewidth=1)                # Border line width

        ax.text(labelXCenterMinutes,
                labelYPosBottomOfText,
                finalLabelTextForGroup, # This is the numerical ID
                ha='center',
                va='bottom', 
                color=groupLabelColor, # Text color (e.g., lime, as set in main config)
                fontsize=groupLabelFontSize,
                weight='bold',
                clip_on=False,
                bbox=bbox_props # Apply the bounding box
                )
        
        # Adjust yMaxFromGroupLabels considering the bbox might affect height.
        # For simplicity, we assume estimatedLineHeightForGroupLabel roughly accounts for it,
        # or that the bbox padding is small enough not to drastically change the overall height calc.
        # If precise bbox height is needed, it's more complex to calculate from text properties.
        textBlockTotalHeightDataUnits = numLinesInLabel * estimatedLineHeightForGroupLabel 
        topOfThisLabelYCoord = labelYPosBottomOfText + textBlockTotalHeightDataUnits
        
        yMaxFromGroupLabels = max(yMaxFromGroupLabels, topOfThisLabelYCoord + 0.1) 
        
    return yMaxFromGroupLabels, group_legend_details # Return legend details
def finalize_plot_styling(fig, ax, timeDf, displayNameMap, yMaxForPlot, 
                          textColor, elementColor, backgroundColor, defaultFontSize,
                          group_legend_details, groupLabelColor, groupBoxEdgeColor): # Added parameters
    """Applies final styling to the plot (labels, ticks, limits, grid, spines, and group legend)."""
    numFunctions = len(timeDf)

    ax.set_xlabel('Time (minutes)', fontsize=defaultFontSize, labelpad=15, color=textColor)
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
            
    xLimUpperMinutes = 1.0 
    if maxEndTimeValSeconds > 0:
        xLimUpperMinutes = (maxEndTimeValSeconds / 60) * 1.05 
    
    ax.set_xlim(0, xLimUpperMinutes)
    
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

    # Add legend for groups
    if group_legend_details:
        legend_elements = []
        for group_id_str, full_group_name in group_legend_details:
            legend_elements.append(
                patches.Patch(facecolor=groupLabelColor, # Color for the patch in legend
                              edgecolor=groupBoxEdgeColor, # Edge color for the patch
                              label=f"{group_id_str}: {full_group_name}")
            )
        
        ax.legend(handles=legend_elements,
                  title="Function Groups",
                  fontsize=defaultFontSize * 0.85,
                  title_fontsize=defaultFontSize * 0.9,
                  loc='upper left', 
                  bbox_to_anchor=(1.02, 1.0), # Position legend outside plot area
                  facecolor=backgroundColor, # Legend box background color
                  edgecolor=elementColor,    # Legend box border color
                  labelcolor=textColor       # Text color for legend items
                 )

    # Adjust layout, especially 'right' to make space for the legend.
    fig.subplots_adjust(left=0.22, right=0.80, top=0.88, bottom=0.08) 

def generate_gantt_chart(config):
    """
    Generates and saves the Gantt chart based on provided function data dictionary
    and saves it to a path relative to reporterDir.
    """
    ## unpack config
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    timeInfo  = config["runtimeInfo"]["timeInfo"]
    # --- Configuration ---
    backgroundColor = 'black'
    textColor = 'yellow'
    elementColor = 'yellow'
    groupLabelColor = 'lime'
    groupBoxEdgeColor = 'lime'
    defaultFontSize = 16
    groupLabelFontSize = 16
    groupLabelMaxCharsPerLine = 25 # Kept in config, passed to draw_group_annotations, but not used for wrapping there
    estimatedLineHeightForGroupLabel = 0.4
    groupBoxPaddingY = 0.10
    groupLabelYOffset = 0.1
    minFigHeight = 9
    figHeightPerFunction = 0.5
    figHeightBasePadding = 2
    # --- End Configuration ---

    # 1. Data Preparation
    data_for_df = []
    for func_name_key, details_dict in timeInfo.items():
        data_for_df.append({
            'functionName': func_name_key,
            'executionTimeSeconds': details_dict['executionTime']
        })
    timeDf = pd.DataFrame(data_for_df)

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
    for func_name_key, details_dict in timeInfo.items():
        csvName = func_name_key
        displayName = details_dict['functionAlias']
        groupName = details_dict['functionGroup']
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
    # Call modified draw_group_annotations, get legend details
    yMaxFromGroupLabels, group_legend_details = draw_group_annotations(ax, timeDf, functionMap, numFunctions, 
                                               groupLabelColor, groupBoxEdgeColor, 
                                               groupBoxPaddingY, groupLabelYOffset, 
                                               groupLabelMaxCharsPerLine, # Param passed but not used for wrapping
                                               estimatedLineHeightForGroupLabel, 
                                               groupLabelFontSize)
    yMaxForPlot = max(yMaxForPlot, yMaxFromGroupLabels)
    
    # Call modified finalize_plot_styling, pass legend details and colors
    finalize_plot_styling(fig, ax, timeDf, displayNameMap, yMaxForPlot, 
                          textColor, elementColor, backgroundColor, defaultFontSize,
                          group_legend_details, groupLabelColor, groupBoxEdgeColor)

    ganttPngPath = p.join(timeImagesDir, "gantt_chart.png")
    # Ensure bbox_inches='tight' to include the legend if it's outside the main axes
    plt.savefig(ganttPngPath, dpi=300, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close(fig)

    relativePath = p.relpath(ganttPngPath, reporterDir)
    return relativePath



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

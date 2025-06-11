import pandas as pd
import os
from os import path as p
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as mpe
import matplotlib.patches as patches
import seaborn as sns
import matplotlib.text as mtext # For creating text artists
from matplotlib.legend_handler import HandlerBase # Base class for custom handlers
import matplotlib.ticker as ticker

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
        "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF",
        "#FF4500", "#32CD32", "#1E90FF", "#FFD700", "#FF69B4", "#00CED1",
        "#FF6347", "#ADFF2F", "#4682B4", "#FFA500", "#DA70D6", "#20B2AA",
        "#DC143C", "#7CFC00"
    ]
    return vibrantColors

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

def draw_task_bars(ax, timeDf, barColors, totalRuntimeSeconds, textColor, backgroundColor, elementColor, defaultFontSize, time_unit_divisor):
    """Draws the Gantt chart bars for each task, scaled by time_unit_divisor."""
    numFunctions = len(timeDf)
    for i, row in timeDf.iterrows():
        yPos = numFunctions - 1 - i
        bar_start_units = row['startTimeSeconds'] / time_unit_divisor
        bar_duration_units = row['executionTimeSeconds'] / time_unit_divisor
        
        ax.broken_barh(
            xranges=[(bar_start_units, bar_duration_units)],
            yrange=(yPos - 0.4, 0.8),
            facecolors=[barColors[i % len(barColors)] if numFunctions > 0 else 'blue'],
            edgecolor=elementColor,
            linewidth=1
        )
        # Display execution time label on significant bars
        if totalRuntimeSeconds > 0 and row['executionTimeSeconds'] > (0.05 * totalRuntimeSeconds):
            timeLabel = format_time_hms(row['executionTimeSeconds'])
            text_x_pos_units = bar_start_units + bar_duration_units / 2
            ax.text(text_x_pos_units, yPos, timeLabel,
                    ha='center', va='center', color=textColor, fontsize=defaultFontSize,
                    path_effects=[mpe.Stroke(linewidth=2, foreground=backgroundColor), mpe.Normal()])

def draw_group_annotations(ax, timeDf, functionMap, numFunctions,
                           groupLabelColor, groupBoxEdgeColor,
                           groupBoxPaddingY, groupLabelYOffset, 
                           groupLabelMaxCharsPerLine, # Unused in current implementation snippet
                           estimatedLineHeightForGroupLabel, # Unused in current implementation snippet
                           groupLabelFontSize,
                           time_unit_divisor, group_label_x_offset):
    """Draws group boxes and numerical labels, scaled by time_unit_divisor."""
    yMaxFromGroupLabels = 0.0
    groupLegendDetails = [] 

    if numFunctions == 0:
        yMaxFromGroupLabels = 1.0 # Default yMax if no functions

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

        rect_x_start_units = boxXStartSeconds / time_unit_divisor
        rect_width_units = (boxXEndSeconds - boxXStartSeconds) / time_unit_divisor
        if rect_width_units < 0: rect_width_units = 0

        rect = patches.Rectangle(
            (rect_x_start_units, rectYStart), rect_width_units, rectHeight,
            linewidth=0.75, edgecolor=groupBoxEdgeColor, facecolor='none', zorder=0.5
        )
        ax.add_patch(rect)

        label_x_pos_units = rect_x_start_units + rect_width_units + group_label_x_offset
        labelYPos = rectYStart + rectHeight / 2

        groupIdToDisplay = str(groupIdx + 1)
        groupLegendDetails.append((groupIdToDisplay, groupName)) 

        bboxProps = dict(boxstyle="square,pad=0.2", facecolor="black", edgecolor="lime", linewidth=1)                
        ax.text(label_x_pos_units, labelYPos, groupIdToDisplay, 
                ha='left', va='center', color=groupLabelColor, fontsize=groupLabelFontSize,
                weight='bold', clip_on=False, bbox=bboxProps)
        
        yMaxFromGroupLabels = max(yMaxFromGroupLabels, boxTopYCoord + groupLabelYOffset)

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
            width / 2.0, height / 2.0, groupIdStr,
            color=self.textColorName, fontsize=fontsize * self.fontSizeFactor, 
            ha="center", va="center", fontweight=self.fontWeight,
            fontfamily=self.fontFamily, transform=trans 
        )
        bboxProps = dict(boxstyle=self.boxStyleStr, facecolor=self.boxFacecolorName,
                         edgecolor=self.boxEdgecolorName, linewidth=self.boxLinewidth)
        textArtist.set_bbox(bboxProps)
        return [textArtist]

def finalize_plot_styling(fig, ax, timeDf, displayNameMap, yMaxForPlot, 
                          textColor, elementColor, backgroundColor, defaultFontSize,
                          groupLegendDetails, groupLabelColor, groupBoxEdgeColor,
                          x_axis_label_text): # Added x_axis_label_text
    """Applies final styling to the plot."""
    numFunctions = len(timeDf)

    ax.set_title('Function Execution Gantt Chart', fontsize=24, pad=20, weight='bold', color=textColor)
    ax.set_xlabel(x_axis_label_text, fontsize=defaultFontSize + 2, color=textColor, labelpad=15, weight='bold') # Set x-axis label

    if numFunctions > 0:
        ax.set_yticks(np.arange(numFunctions))
        originalYTickLabels = timeDf['functionName'].iloc[::-1].tolist() 
        displayYTickLabels = [displayNameMap.get(name, name) for name in originalYTickLabels]
        ax.set_yticklabels(displayYTickLabels, fontsize=defaultFontSize)
    else:
        ax.set_yticks([])

    currentYLimBottom, _ = ax.get_ylim() 
    ax.set_ylim(currentYLimBottom, yMaxForPlot)

    ax.grid(True, axis='x', linestyle=':', alpha=0.6, color=elementColor)
    ax.grid(False, axis='y') 
    ax.tick_params(axis='x', colors=textColor, labelsize=defaultFontSize)
    ax.tick_params(axis='y', colors=textColor, labelsize=defaultFontSize) 
    for spine_pos in ['bottom', 'left']:
        ax.spines[spine_pos].set_color(elementColor)
    for spine_pos in ['top', 'right']:
        ax.spines[spine_pos].set_color(backgroundColor)
    sns.despine(fig=fig, ax=ax, top=True, right=True, left=False, bottom=False) # Pass fig and ax to despine

    if groupLegendDetails: 
        legendHandles = []
        legendLabels = []
        
        styledIdKeyHandler = HandlerStyledIdBox(
            textColorName=groupLabelColor, boxFacecolorName='black', 
            boxEdgecolorName=groupLabelColor, # Edge of ID box matches label color
            fontSizeFactor=0.85, boxStyleStr='square,pad=0.3', boxLinewidth=1
        )
        handlerMapForLegend = {}
        for groupIdStr, fullGroupName in groupLegendDetails:
            legendHandles.append(groupIdStr)
            legendLabels.append(f"{fullGroupName}") 
            handlerMapForLegend[groupIdStr] = styledIdKeyHandler

        legend_ncol = 1
        if len(legendHandles) > 5: # Example: use multiple columns for many groups
            legend_ncol = max(1, len(legendHandles) // 5 + (1 if len(legendHandles) % 5 > 0 else 0) )


        ax.legend(
            legendHandles, legendLabels,   
            title="Groups:", # Added legend title
            fontsize=defaultFontSize * 0.9, title_fontsize=defaultFontSize * 0.9,
            loc='lower left', bbox_to_anchor=(0.0, 0.0), # Original position
            facecolor=backgroundColor, edgecolor=groupLabelColor, # Legend box edge color (lime)
            labelcolor=textColor, handler_map=handlerMapForLegend,
            framealpha=0.8, ncol=legend_ncol
        )
    
    # Adjust bottom more if legend is present to avoid overlap with x-label
    bottom_padding = 0.20 if groupLegendDetails else 0.15
    fig.subplots_adjust(left=0.22, right=0.95, top=0.90, bottom=bottom_padding)


# Main Visualization Orchestration Functions
def generate_gantt_chart(config):
    """
    Generates and saves the Gantt chart. X-axis is in hours if total runtime > 5 hours, else minutes.
    """
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    timeInfo  = config["runtimeInfo"]["timeInfo"]
    
    # --- Style Configuration ---
    backgroundColor = 'black'
    textColor = 'yellow'
    elementColor = 'yellow'
    groupLabelColor = 'lime'
    groupBoxEdgeColor = 'lime' # For the rectangle around grouped tasks
    defaultFontSize = 16
    groupLabelFontSize = 16 # For the numeric ID next to group boxes
    
    # Unused parameters from original snippet, kept for potential future use or if missed elsewhere
    groupLabelMaxCharsPerLine = 25 
    estimatedLineHeightForGroupLabel = 0.4
    
    groupBoxPaddingY = 0.10 # Vertical padding for group boxes
    groupLabelYOffset = 0.1 # Vertical offset for calculating yMax based on group labels
    minFigHeight = 9.0
    figHeightPerFunction = 0.5
    figHeightBasePadding = 2.0

    # --- Data Preparation ---
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
        # Assuming the input order from timeInfo defines the plotting order of functions (top to bottom)
        # If a specific sort order is needed, apply it here before calculating cumulative sums.
        # E.g., timeDf = timeDf.sort_values(by='startTimeThemselves', ascending=True).reset_index(drop=True)
        timeDf['startTimeSeconds'] = timeDf['executionTimeSeconds'].cumsum().shift(1).fillna(0)
        timeDf['endTimeSeconds'] = timeDf['startTimeSeconds'] + timeDf['executionTimeSeconds']
        totalRuntimeSeconds = timeDf['executionTimeSeconds'].sum()
        if pd.isna(totalRuntimeSeconds): totalRuntimeSeconds = 0.0

    displayNameMap = {} # Maps internal functionName to display alias
    functionMap = {}    # Maps groupName to {functionName: functionAlias}
    for funcNameKey, detailsDict in timeInfo.items():
        csvName = funcNameKey 
        displayName = detailsDict['functionAlias']
        groupName = detailsDict['functionGroup']
        displayNameMap[csvName] = displayName
        if groupName not in functionMap:
            functionMap[groupName] = {}
        functionMap[groupName][csvName] = displayName

    # --- Determine X-axis Time Unit and Label Offset ---
    RUNTIME_THRESHOLD_SECONDS = 5 * 3600 # 5 hours
    if totalRuntimeSeconds > RUNTIME_THRESHOLD_SECONDS:
        time_unit_divisor = 3600.0  # Hours
        x_axis_label_text = "Time (Hours)"
    else:
        time_unit_divisor = 60.0   # Minutes
        x_axis_label_text = "Time (Minutes)"

    max_x_val_on_axis = totalRuntimeSeconds / time_unit_divisor if totalRuntimeSeconds > 0 else 0
    
    # Relative offset for group ID labels from their boxes
    if max_x_val_on_axis > 1e-9: # Check for non-zero scaled runtime
        group_label_x_offset = 0.005 * max_x_val_on_axis 
    else: # Fallback for zero or very small total runtime
        group_label_x_offset = 0.1 if time_unit_divisor == 60.0 else 0.01 # 0.1 min or 0.01 hr

    # --- Plotting ---
    timeImagesDir = p.join(reporterDir, "Images", "time_data")
    os.makedirs(timeImagesDir, exist_ok=True)

    initialize_matplotlib_settings(defaultFontSize)
    numFunctions = len(timeDf)
    fig, ax = create_figure_and_axes(numFunctions, backgroundColor, minFigHeight, figHeightPerFunction, figHeightBasePadding)
    
    # Use a vibrant color palette, ensuring at least one color if numFunctions is 0 (though no bars would be drawn)
    barColors = sns.color_palette('viridis', n_colors=max(1, numFunctions))

    draw_task_bars(ax, timeDf, barColors, totalRuntimeSeconds, textColor, backgroundColor, elementColor, defaultFontSize, time_unit_divisor)
    
    yMaxForPlot = numFunctions if numFunctions > 0 else 1.0 # Base y-limit
    yMaxFromGroupLabels, groupLegendDetails = draw_group_annotations(
        ax, timeDf, functionMap, numFunctions, 
        groupLabelColor, groupBoxEdgeColor, 
        groupBoxPaddingY, groupLabelYOffset, 
        groupLabelMaxCharsPerLine, estimatedLineHeightForGroupLabel, 
        groupLabelFontSize,
        time_unit_divisor, group_label_x_offset
    )
    yMaxForPlot = max(yMaxForPlot, yMaxFromGroupLabels) # Adjust y-limit for group annotations
    
    finalize_plot_styling(
        fig, ax, timeDf, displayNameMap, yMaxForPlot, 
        textColor, elementColor, backgroundColor, defaultFontSize,
        groupLegendDetails, groupLabelColor, groupBoxEdgeColor, # Pass groupBoxEdgeColor for consistency
        x_axis_label_text
    )

    ganttPngPath = p.join(timeImagesDir, "gantt_chart.png")
    plt.savefig(ganttPngPath, dpi=300, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close(fig)

    relativePath = p.relpath(ganttPngPath, reporterDir)
    return relativePath
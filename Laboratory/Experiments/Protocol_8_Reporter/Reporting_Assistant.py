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
# Note: matplotlib.pyplot is imported as plt earlier in the file, no need to re-import
# import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import matplotlib.patheffects as mpe
import numpy as np
import seaborn as sns
# import textwrap # Already imported
import os
import os.path as p
from typing import List, Dict, Tuple, Union, Any # Added for type hinting

# Placeholder classes (extend if needed)
class FilePath:
    pass
class DirectoryPath:
    pass

def format_time_hms(seconds: Union[float, int]) -> str:
    """Converts seconds to HH:MM:SS string format."""
    seconds = int(round(seconds))
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    secs = seconds % 60
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"

def initialize_matplotlib_settings(default_font_size: int) -> None:
    """Sets global Matplotlib rcParams and verifies font."""
    plt.rcParams['font.family'] = 'monospace'
    plt.rcParams['font.monospace'] = ['Consolas', 'DejaVu Sans Mono', 'Courier New'] # type: ignore
    plt.rcParams['font.size'] = default_font_size

# load_and_prepare_data is removed as its functionality for the new input
# format is integrated into generate_gantt_chart.

def create_figure_and_axes(num_functions: int, background_color: str, 
                           min_fig_height: float, fig_height_per_function: float, 
                           fig_height_base_padding: float) -> Tuple[plt.Figure, plt.Axes]:
    """Creates the Matplotlib figure and axes with basic styling."""
    figHeight = max(min_fig_height, num_functions * fig_height_per_function + fig_height_base_padding)
    fig, ax = plt.subplots(figsize=(24, figHeight))
    fig.patch.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    return fig, ax

def draw_task_bars(ax: plt.Axes, time_df: pd.DataFrame, bar_colors: List[str], 
                   total_runtime_seconds: float, text_color: str, background_color: str, 
                   element_color: str, default_font_size: int) -> None:
    """Draws the Gantt chart bars for each task."""
    numFunctions = len(time_df)
    for i, row in time_df.iterrows():
        yPos = numFunctions - 1 - i
        barStartMinutes = row['startTimeSeconds'] / 60
        barDurationMinutes = row['executionTimeSeconds'] / 60
        
        ax.broken_barh(
            xranges=[(barStartMinutes, barDurationMinutes)],
            yrange=(yPos - 0.4, 0.8),
            facecolors=[bar_colors[i % len(bar_colors)] if numFunctions > 0 else 'blue'], # type: ignore
            edgecolor=element_color,
            linewidth=0.75
        )
        if total_runtime_seconds > 0 and row['executionTimeSeconds'] > (0.05 * total_runtime_seconds):
            timeLabel = format_time_hms(row['executionTimeSeconds'])
            textXPosMinutes = barStartMinutes + barDurationMinutes / 2
            ax.text(textXPosMinutes,
                    yPos, timeLabel,
                    ha='center', va='center', color=text_color, fontsize=default_font_size,
                    path_effects=[mpe.Stroke(linewidth=2, foreground=background_color), mpe.Normal()])

def draw_group_annotations(ax: plt.Axes, time_df: pd.DataFrame, function_map: Dict[str, Dict[str, str]], 
                           num_functions: int, group_label_color: str, group_box_edge_color: str, 
                           group_box_padding_y: float, group_label_y_offset: float, 
                           group_label_max_chars_per_line: int, estimated_line_height_for_group_label: float, 
                           group_label_font_size: int) -> float:
    """Draws group boxes and labels on the Gantt chart."""
    yMaxFromGroupLabels = 0.0
    if num_functions == 0: 
        yMaxFromGroupLabels = 1.0 # Default if no functions, to ensure plot has some height.

    for groupName, functionsInGroupMap in function_map.items():
        csvNamesInGroup = list(functionsInGroupMap.keys())
        groupMemberRows = time_df[time_df['functionName'].isin(csvNamesInGroup)]
        if groupMemberRows.empty:
            continue

        memberIndices = groupMemberRows.index.tolist()
        if not memberIndices: # Should not happen if groupMemberRows is not empty
            continue

        currentMinIndexInDf = min(memberIndices)
        currentMaxIndexInDf = max(memberIndices)
        
        yPosOfTopBarInGroup = num_functions - 1 - currentMinIndexInDf
        yPosOfBottomBarInGroup = num_functions - 1 - currentMaxIndexInDf

        boxBottomYCoord = yPosOfBottomBarInGroup - 0.4 - group_box_padding_y
        boxTopYCoord = yPosOfTopBarInGroup + 0.4 + group_box_padding_y
        
        rectYStart = boxBottomYCoord
        rectHeight = boxTopYCoord - rectYStart
        if rectHeight <= 0: continue

        boxXStartSeconds = groupMemberRows['startTimeSeconds'].min()
        boxXEndSeconds = groupMemberRows['endTimeSeconds'].max()
        
        rectXStartMinutes = boxXStartSeconds / 60
        rectWidthMinutes = (boxXEndSeconds - boxXStartSeconds) / 60
        if rectWidthMinutes < 0: rectWidthMinutes = 0 # Ensure non-negative width

        rect = patches.Rectangle(
            (rectXStartMinutes, rectYStart),
            rectWidthMinutes,
            rectHeight,
            linewidth=1.5,
            edgecolor=group_box_edge_color,
            facecolor='none',
            zorder=0.5 # Draw behind bars slightly if overlap, though ideally not
        )
        ax.add_patch(rect)

        labelXCenterMinutes = rectXStartMinutes + rectWidthMinutes / 2
        labelYPosBottomOfText = boxTopYCoord + group_label_y_offset
        
        finalLabelTextForGroup = groupName
        if len(groupName) > group_label_max_chars_per_line:
            finalLabelTextForGroup = textwrap.fill(groupName, width=group_label_max_chars_per_line)
        
        numLinesInLabel = finalLabelTextForGroup.count('\n') + 1
        
        ax.text(labelXCenterMinutes,
                labelYPosBottomOfText,
                finalLabelTextForGroup,
                ha='center',
                va='bottom', # Anchor text at its bottom
                color=group_label_color,
                fontsize=group_label_font_size,
                weight='bold',
                clip_on=False # Allow labels to go outside plot area if necessary before final adjustment
                )
        
        textBlockTotalHeightDataUnits = numLinesInLabel * estimated_line_height_for_group_label # Data units
        topOfThisLabelYCoord = labelYPosBottomOfText + textBlockTotalHeightDataUnits
        
        yMaxFromGroupLabels = max(yMaxFromGroupLabels, topOfThisLabelYCoord + 0.1) # Add small padding
        
    return yMaxFromGroupLabels


def finalize_plot_styling(fig: plt.Figure, ax: plt.Axes, time_df: pd.DataFrame, 
                          display_name_map: Dict[str, str], y_max_for_plot: float, 
                          text_color: str, element_color: str, background_color: str, 
                          default_font_size: int) -> None:
    """Applies final styling to the plot (labels, ticks, limits, grid, spines)."""
    numFunctions = len(time_df)

    ax.set_xlabel('Time (minutes)', fontsize=default_font_size, labelpad=15, color=text_color)
    ax.set_title('Function Execution Gantt Chart', fontsize=24, pad=20, weight='bold', color=text_color)

    if numFunctions > 0:
        ax.set_yticks(np.arange(numFunctions))
        # Order of y-tick labels should correspond to order in timeDf, plotted from top to bottom
        originalYTickLabels = time_df['functionName'].iloc[::-1].tolist() # Reversed for set_yticklabels
        displayYTickLabels = [display_name_map.get(name, name) for name in originalYTickLabels]
        ax.set_yticklabels(displayYTickLabels, fontsize=default_font_size)
    else:
        ax.set_yticks([])

    maxEndTimeValSeconds = 0.0
    if not time_df.empty:
        maxEndTimeValSeconds = time_df['endTimeSeconds'].max()
        if pd.isna(maxEndTimeValSeconds): # handle case where max is NaN
            maxEndTimeValSeconds = 0.0
            
    xLimUpperMinutes = 1.0 # Default if no tasks or all zero duration
    if maxEndTimeValSeconds > 0:
        xLimUpperMinutes = (maxEndTimeValSeconds / 60) * 1.05 # Add 5% padding
    
    ax.set_xlim(0, xLimUpperMinutes)
    
    currentYLimBottom, _ = ax.get_ylim() # Preserve default bottom limit (e.g., -0.5)
    ax.set_ylim(currentYLimBottom, y_max_for_plot)

    ax.grid(True, axis='x', linestyle=':', alpha=0.6, color=element_color)
    ax.grid(False, axis='y') # No horizontal grid lines for y-axis
    ax.tick_params(axis='x', colors=text_color, labelsize=default_font_size)
    ax.tick_params(axis='y', colors=text_color, labelsize=default_font_size) # Color for y-tick labels
    ax.spines['bottom'].set_color(element_color)
    ax.spines['left'].set_color(element_color)
    ax.spines['top'].set_color(background_color) # Effectively hide top spine
    ax.spines['right'].set_color(background_color) # Effectively hide right spine
    sns.despine(top=True, right=True, left=False, bottom=False) # Keep left and bottom spines visible

    fig.subplots_adjust(left=0.22, right=0.96, top=0.88, bottom=0.08) # Adjust layout

def generate_gantt_chart(config: dict) -> FilePath:
    """
    Generates and saves the Gantt chart based on provided function data dictionary
    and saves it to a path relative to reporterDir.
    """
    ## unpack config
    reporterDir: DirectoryPath = config["runtimeInfo"]["madeByReporting"]["reporterDir"] # type: ignore
    timeInfo: Dict[str, Dict[str, Any]]  = config["runtimeInfo"]["timeInfo"]
    # --- Configuration ---
    backgroundColor = 'black'
    textColor = 'yellow'
    elementColor = 'yellow'
    groupLabelColor = 'lime'
    groupBoxEdgeColor = '#777777'
    defaultFontSize = 16
    groupLabelFontSize = 16
    groupLabelMaxCharsPerLine = 25
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
        # Initialize empty DataFrame with expected columns if input is empty
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

    # Directory setup
    timeImagesDir = p.join(reporterDir, "Images", "time_data")
    os.makedirs(timeImagesDir, exist_ok=True)

    # 2. Plotting
    initialize_matplotlib_settings(default_font_size=defaultFontSize)
    numFunctions = len(timeDf)
    fig, ax = create_figure_and_axes(num_functions=numFunctions, background_color=backgroundColor, 
                                     min_fig_height=minFigHeight, fig_height_per_function=figHeightPerFunction, 
                                     fig_height_base_padding=figHeightBasePadding)
    
    barColors = sns.color_palette('viridis', n_colors=max(1, numFunctions))

    draw_task_bars(ax=ax, time_df=timeDf, bar_colors=barColors, total_runtime_seconds=totalRuntimeSeconds, 
                   text_color=textColor, background_color=backgroundColor, element_color=elementColor, 
                   default_font_size=defaultFontSize)
    
    yMaxForPlot = float(numFunctions) if numFunctions > 0 else 1.0
    yMaxFromGroupLabels = draw_group_annotations(ax=ax, time_df=timeDf, function_map=functionMap, num_functions=numFunctions, 
                                               group_label_color=groupLabelColor, group_box_edge_color=groupBoxEdgeColor, 
                                               group_box_padding_y=groupBoxPaddingY, group_label_y_offset=groupLabelYOffset, 
                                               group_label_max_chars_per_line=groupLabelMaxCharsPerLine, 
                                               estimated_line_height_for_group_label=estimatedLineHeightForGroupLabel, 
                                               group_label_font_size=groupLabelFontSize)
    yMaxForPlot = max(yMaxForPlot, yMaxFromGroupLabels)
    
    finalize_plot_styling(fig=fig, ax=ax, time_df=timeDf, display_name_map=displayNameMap, y_max_for_plot=yMaxForPlot, 
                          text_color=textColor, element_color=elementColor, background_color=backgroundColor, 
                          default_font_size=defaultFontSize)

    ganttPngPath: FilePath = p.join(timeImagesDir, "gantt_chart.png") # type: ignore
    plt.savefig(ganttPngPath, dpi=300, facecolor=fig.get_facecolor(), bbox_inches='tight') # type: ignore
    plt.close(fig)

    relativePath: FilePath = p.relpath(ganttPngPath, reporterDir) # type: ignore
    return relativePath



def make_charge_visualisation(config: dict, out_dir: DirectoryPath) -> Tuple[FilePath, float, float]:
    ## unpack config
    cappedPdb: FilePath = config["runtimeInfo"]["madeByCapping"]["cappedPdb"] # type: ignore
    chargesCsv: FilePath = config["runtimeInfo"]["madeByCharges"]["chargesCsv"] # type: ignore
    reporterDir: DirectoryPath = config["runtimeInfo"]["madeByReporting"]["reporterDir"] # type: ignore
    cappedDf = pdbUtils.pdb2df(cappedPdb)
    chargesDf = pd.read_csv(chargesCsv, index_col="Unnamed: 0")

    chargeGradientColors = make_charge_gradient_colors()

    cappedDf["BETAFACTOR"] = chargesDf["Charge"]

    annotatedPdb: FilePath = p.join(out_dir, "annotated_charged.pdb") # type: ignore
    pdbUtils.df2pdb(cappedDf, annotatedPdb)

    minCharge: float = chargesDf["Charge"].min()
    maxCharge: float = chargesDf["Charge"].max()

    ## load PDB into a string
    with open(annotatedPdb, 'r') as f:
        pdbBlock = "".join(f.readlines())


    # Create 3D visualization with py3Dmol
    view = py3Dmol.view(width=600, height=600) # Increased size slightly for better viewing
    view.addModel(pdbBlock, 'pdb', {"keepH": True}) # type: ignore
    view.setStyle({ # type: ignore
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
    view.setBackgroundColor("black") # type: ignore

    for index, row in chargesDf.iterrows():
        view.addLabel(row["ATOM_NAME"], # type: ignore
                        {'fontSize': 10, 'fontColor': 'yellow', 
                        'showBackground': True, 'backgroundOpacity': 1, 'backgroundColor': 'black',
                        'alignment': 'topCenter'},
                        {'serial': index+1}) 
        view.addLabel(f"{round(row['Charge'],2):.2f}", # type: ignore
                        {'fontSize': 10, 'fontColor': 'yellow', 
                        'showBackground': True, 'backgroundOpacity': 1, 'backgroundColor': 'black',
                        'alignment': 'bottomCenter'},
                        {'serial': index+1}) 
    view.zoomTo() # type: ignore

    # Save the visualization as HTML
    chargeHtml: FilePath = p.join(out_dir, "charges.html") # type: ignore
    with open(chargeHtml, 'w') as f:
        f.write(view._make_html()) # type: ignore

    os.remove(annotatedPdb)

    relativePath: FilePath = p.relpath(chargeHtml, reporterDir) # type: ignore

    return relativePath, minCharge, maxCharge


def make_charge_color_bar(min_charge: float, max_charge: float, out_dir: DirectoryPath, config: dict) -> FilePath:
    reporterDir: DirectoryPath = config["runtimeInfo"]["madeByReporting"]["reporterDir"] # type: ignore
    chargeGradientColors = make_charge_gradient_colors() 

    # Create figure with black background
    fig, ax = plt.subplots(figsize=(6, 1), facecolor='black') # Figure size in inches
    fig.subplots_adjust(bottom=0.5) # Adjust subplot to leave space for label below color bar

    cmap = mcolors.ListedColormap(chargeGradientColors)
    norm = mcolors.Normalize(vmin=min_charge, vmax=max_charge)

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

    outputFilePath: FilePath = os.path.join(out_dir, "charge_color_bar.png") # type: ignore
    
    # Save the figure with a specific DPI to control pixel height
    plt.savefig(outputFilePath, facecolor='black', edgecolor='none') # type: ignore
    plt.close(fig) # Close the figure to free memory

    relativePath: FilePath = p.relpath(outputFilePath, reporterDir) # type: ignore
    return relativePath


def  make_conformer_visualisations(config: dict, out_dir: DirectoryPath) -> List[FilePath]:
    ## unpack config
    conformerXyzs: List[FilePath] = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"] # type: ignore
    reporterDir: DirectoryPath =     config["runtimeInfo"]["madeByReporting"]["reporterDir"] # type: ignore

    vibrantColors = make_vibrant_colors()

    conformerHtmls: List[FilePath] = []
    for idx, xyzFile_item in enumerate(conformerXyzs): # Renamed loop var
        name = xyzFile_item.split("/")[-1].split(".")[0]
        xyzBlock = xyz_to_string(xyz_file=xyzFile_item)
        # Create 3D visualization with py3Dmol
        view = py3Dmol.view(width=200, height=200) 
        view.addModel(xyzBlock, 'xyz', {"keepH": True}) # type: ignore
        # Use modulo to cycle through colors
        color = vibrantColors[idx % len(vibrantColors)]
        view.setStyle({"elem": "C"},{ # type: ignore
            'stick': {
                'radius': 0.2, "color": color
            },
            'sphere': {
                'scale': 0.3, "color": color
            }
        })
        view.setStyle({"elem": "C", "invert": True},{ # type: ignore
            'stick': {
                'radius': 0.2
            },
            'sphere': {
                'scale': 0.3
            }
        })

        view.setBackgroundColor("black") # type: ignore

        # Save the visualization as HTML
        conformerHtml: FilePath = p.join(out_dir, f"{name}.html") # type: ignore
        with open(conformerHtml, 'w') as f:
            f.write(view._make_html()) # type: ignore

        relativePath: FilePath = p.relpath(conformerHtml, reporterDir) # type: ignore

        conformerHtmls.append(relativePath)

    return conformerHtmls


def make_highlighted_torsion_visualisations(config: dict, out_dir: DirectoryPath) -> Dict[str, FilePath]:
    ## unpack config
    cappedPdb: FilePath = config["runtimeInfo"]["madeByCapping"]["cappedPdb"] # type: ignore
    uniqueRotatableDihedrals: Dict[str, Dict[str, List[Any]]] = config["runtimeInfo"]["madeByTwisting"]["rotatableDihedrals"]
    reporterDir: DirectoryPath =     config["runtimeInfo"]["madeByReporting"]["reporterDir"] # type: ignore

    pdbBlock = pdb_to_string(pdb_file=cappedPdb)

    htmlFiles: Dict[str, FilePath] = {}
    for torsionTag, dihedralData in uniqueRotatableDihedrals.items():
        view = py3Dmol.view(width=300, height=300)
        view.addModel(pdbBlock, "pdb",  {"keepH": True}) # type: ignore

        ## simple whiteCarbons balls and sticks
        view.setStyle({},{ # type: ignore
            'stick': {
                'radius': 0.2
            },
            'sphere': {
                'scale': 0.3
            }
        })
        view.setStyle({"atom": dihedralData["ATOM_NAMES"]}, { # type: ignore
            'sphere': {
                'scale': 0.3,
                'color': "#00FF00"
            },
                    'stick': {
                'radius': 0.2
            }
        })
        view.setBackgroundColor("black") # type: ignore

        # Save the visualization as HTML
        htmlFile: FilePath = p.join(out_dir, f"{torsionTag}.html") # type: ignore
        with open(htmlFile, 'w') as f:
            f.write(view._make_html()) # type: ignore

        relativePath: FilePath = p.relpath(htmlFile, reporterDir) # type: ignore
        htmlFiles[torsionTag] = relativePath

    return htmlFiles

def make_vibrant_colors() -> List[str]:
    """Returns a list of vibrant hex color codes."""
    vibrantColors = [ # Renamed from vibrantColours
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
    return vibrantColors


def make_charge_gradient_colors() -> List[str]: # Renamed from make_charge_gradient_colours
    """Returns a list of hex color codes for a charge gradient."""
    chargeGradientColors = ['#FF00FF', '#FF66FF', '#FFCCFF', '#FFF5FF', '#FFFFFF', '#F5FFF5', '#CCFFCC', '#33FF33', '#00FF00']
    return chargeGradientColors

def pdb_to_string(pdb_file: FilePath) -> str: # Renamed argument
    """Reads a PDB file and returns its content as a string."""
    with open(pdb_file, "r") as f: # type: ignore
        pdbBlock = "".join(f.readlines())
    return pdbBlock


def xyz_to_string(xyz_file: FilePath) -> str: # Renamed function and argument
    """Reads an XYZ file and returns its content as a string."""
    with open(xyz_file, 'r') as f: # type: ignore
        xyzBlock = "".join(f.readlines())
    return xyzBlock

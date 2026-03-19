import pandas as pd
import py3Dmol
from pdbUtils import pdbUtils # Your custom PDB parsing utility
from shutil import copy
import os
from os import path as p
from matplotlib import pyplot as plt
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
    vibrantColors =[
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
    chargeGradientColors =['#FF00FF', '#FF66FF', '#FFCCFF', '#FFF5FF', '#FFFFFF', '#F5FFF5', '#CCFFCC', '#33FF33', '#00FF00']
    return chargeGradientColors

def compute_physics_layout(atoms_df, bonds_df, iterations=300, lr=0.05):
    """Simulates springs to flatten molecule while maintaining bonds and angles."""
    import networkx as nx
    import numpy as np
    
    coords = atoms_df[['x', 'y', 'z']].values
    atom_ids = atoms_df['atom_id'].values
    idx_map = {atom_id: i for i, atom_id in enumerate(atom_ids)}

    # Build graph to find relationships
    G = nx.Graph()
    G.add_nodes_from(range(len(atom_ids)))
    for _, row in bonds_df.iterrows():
        if row['atom1_id'] in idx_map and row['atom2_id'] in idx_map:
            G.add_edge(idx_map[row['atom1_id']], idx_map[row['atom2_id']])

    springs =[]
    spring_pairs = set()

    # A. Bond Springs (1-2 pairs) -> Strict Hooke's Law
    for u, v in G.edges():
        d3d = np.linalg.norm(coords[u] - coords[v])
        springs.append((u, v, d3d, 2.0))  # Force constant k=2.0
        spring_pairs.update([(u, v), (v, u)])

    # B. Angle Springs (1-3 pairs) -> Locks geometry (ChemDraw style)
    for node in G.nodes():
        neighbors = list(G.neighbors(node))
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                u, v = neighbors[i], neighbors[j]
                if (u, v) not in spring_pairs:
                    d3d = np.linalg.norm(coords[u] - coords[v])
                    springs.append((u, v, d3d, 1.0)) # Force constant k=1.0
                    spring_pairs.update([(u, v), (v, u)])

    # C. Initialize 2D coordinates using PCA (Find flattest projection to start)
    coords_centered = coords - np.mean(coords, axis=0)
    cov = np.cov(coords_centered.T)
    evals, evecs = np.linalg.eigh(cov)
    pos_2d = coords_centered @ evecs[:, -2:] # Project onto 2 largest principal components
    pos_2d += np.random.normal(scale=0.05, size=pos_2d.shape) # Jitter to prevent locking

    # D. Physics Simulation (Iterative Force Application)
    for _ in range(iterations):
        forces = np.zeros_like(pos_2d)
        
        # Apply Spring Forces
        for u, v, d0, k in springs:
            delta = pos_2d[u] - pos_2d[v]
            d = np.linalg.norm(delta)
            if d > 1e-4:
                force_mag = k * (d - d0)
                force_vec = force_mag * (delta / d)
                forces[u] -= force_vec  # Pull towards ideal
                forces[v] += force_vec
                
        # Apply Repulsive Forces (for non-bonded atoms collapsing during flattening)
        for u in range(len(atom_ids)):
            for v in range(u + 1, len(atom_ids)):
                if (u, v) not in spring_pairs:
                    delta = pos_2d[u] - pos_2d[v]
                    d = np.linalg.norm(delta)
                    if 1e-4 < d < 1.5:  # Only repel if they get suspiciously close
                        force_mag = 0.5 * (1.5 - d) / (d**2)
                        force_vec = force_mag * (delta / d)
                        forces[u] += force_vec  # Push apart
                        forces[v] -= force_vec

        # Update coordinates (Gradient Descent)
        pos_2d += lr * forces

    # E. Transpose to Landscape Orientation
    x_span = np.max(pos_2d[:, 0]) - np.min(pos_2d[:, 0])
    y_span = np.max(pos_2d[:, 1]) - np.min(pos_2d[:, 1])
    
    # If the molecule is taller than it is wide, swap X and Y axes
    if y_span > x_span:
        pos_2d = pos_2d[:, [1, 0]]  # Swaps column 0 and column 1

    return {atom_ids[i]: (pos_2d[i, 0], pos_2d[i, 1]) for i in range(len(atom_ids))}

# Gantt Chart Helper Functions
def initialize_matplotlib_settings(defaultFontSize):
    """Sets global Matplotlib rcParams and verifies font."""
    plt.rcParams['font.family'] = 'monospace'
    plt.rcParams['font.monospace'] =['Consolas', 'DejaVu Sans Mono', 'Courier New']
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
    groupLegendDetails =[] 

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
        displayYTickLabels =[displayNameMap.get(name, name) for name in originalYTickLabels]
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
        legendLabels =[]
        
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

    dataForDf =[]
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

    conformerHtmls =[]
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
    """
    Create Plotly-based 2D visualizations highlighting torsion atoms.
    Uses a custom physics-based spring layout to enforce exact bond lengths 
    and angles, resulting in a ChemDraw-style 2D projection.
    Automatically rotates the coords to ensure a landscape orientation.
    """
    import plotly.graph_objects as go
    import networkx as nx
    import pandas as pd
    import numpy as np
    import os.path as p

    # --- 1. Helper Parsers ---
    def parse_mol2_file(mol2_file):
        with open(mol2_file, 'r') as f:
            lines = f.readlines()

        atom_start = lines.index('@<TRIPOS>ATOM\n') + 1
        bond_start = lines.index('@<TRIPOS>BOND\n') + 1

        atom_lines =[]
        for i in range(atom_start, len(lines)):
            if lines[i].startswith('@<TRIPOS>'): break
            atom_lines.append(lines[i])

        bond_lines =[]
        for i in range(bond_start, len(lines)):
            if lines[i].startswith('@<TRIPOS>'): break
            bond_lines.append(lines[i])

        atoms_df = pd.DataFrame([line.split() for line in atom_lines], 
                                columns=['atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id', 'subst_name', 'charge'])
        for col in['atom_id', 'x', 'y', 'z']:
            atoms_df[col] = pd.to_numeric(atoms_df[col])

        bonds_df = pd.DataFrame([line.split()[:3] for line in bond_lines if line.strip()], 
                                columns=['bond_id', 'atom1_id', 'atom2_id'])
        bonds_df[['bond_id', 'atom1_id', 'atom2_id']] = bonds_df[['bond_id', 'atom1_id', 'atom2_id']].astype(int)

        return atoms_df, bonds_df

    # --- 2. Main Workflow ---
    cappedMol2 = config["runtimeInfo"]["madeByAssembly"]["cappedMol2"]
    uniqueRotatableDihedrals = config["runtimeInfo"]["madeByTwisting"]["torsionsToScan"]
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    
    atoms_df, bonds_df = parse_mol2_file(cappedMol2)
    
    # Run the custom physics simulation to fetch beautifully scaled 2D coords
    pos_2d = compute_physics_layout(atoms_df, bonds_df)
    
    atom_names = atoms_df['atom_name'].values
    
    def get_element(atom_name):
        name = atom_name.strip()
        if len(name) >= 2 and name[1].islower(): return name[:2]
        return name[0]
    
    atom_elements =[get_element(name) for name in atom_names]
    element_colors = {'C': '#808080', 'H': '#FFFFFF', 'O': '#FF0000', 'N': '#0000FF', 'S': '#FFFF00'}
    
    htmlFiles = {}
    for torsionTag, dihedralData in uniqueRotatableDihedrals.items():
        highlighted_indices = dihedralData["ATOM_INDEXES"] 
        
        edge_x, edge_y = [],[]
        for _, bond in bonds_df.iterrows():
            atom1_id, atom2_id = bond['atom1_id'], bond['atom2_id']
            if atom1_id in pos_2d and atom2_id in pos_2d:
                x0, y0 = pos_2d[atom1_id]
                x1, y1 = pos_2d[atom2_id]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])
        
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y, mode='lines',
            line=dict(color='#666666', width=3),
            hoverinfo='none', showlegend=False
        )
        
        node_x =[pos_2d[a_id][0] for a_id in atoms_df['atom_id']]
        node_y = [pos_2d[a_id][1] for a_id in atoms_df['atom_id']]
        node_colors, node_sizes, hover_text, text_labels = [],[], [],[]
        
        for i in range(len(atoms_df)):
            element = atom_elements[i]
            node_colors.append(element_colors.get(element, '#FFA500'))
            node_sizes.append(8 if element == 'H' else 14)
            
            if i in highlighted_indices:
                hover_text.append(f"<b>{atom_names[i]}</b> (Torsion Atom)")
                text_labels.append(atom_names[i])
            else:
                hover_text.append(f"{atom_names[i]} ({element})")
                text_labels.append('')
        
        highlight_x =[node_x[i] for i in highlighted_indices if i < len(node_x)]
        highlight_y = [node_y[i] for i in highlighted_indices if i < len(node_y)]
        
        highlight_trace = go.Scatter(
            x=highlight_x, y=highlight_y, mode='markers',
            marker=dict(size=28, color='#00FF00', opacity=0.3, line=dict(color='#00FF00', width=2)),
            hoverinfo='skip', showlegend=False
        )
        
        node_trace = go.Scatter(
            x=node_x, y=node_y, mode='markers+text',
            marker=dict(size=node_sizes, color=node_colors, line=dict(color='#000000', width=1.5)),
            text=text_labels, textposition='top center',
            textfont=dict(size=10, color='#FFFF00', family='Arial Black'),
            hovertext=hover_text, hoverinfo='text', showlegend=False
        )
        
        fig = go.Figure(data=[edge_trace, highlight_trace, node_trace])
        
        fig.update_layout(
            title=dict(text=f'Torsion: {torsionTag}', font=dict(color='white', size=16)),
            xaxis=dict(visible=False), yaxis=dict(visible=False, scaleanchor='x', scaleratio=1),
            paper_bgcolor='black', plot_bgcolor='black', margin=dict(l=0, r=0, t=40, b=0),
            hovermode='closest', autosize=True, height=330
        )
        
        htmlFile = p.join(outDir, f"{torsionTag}.html")
        fig.write_html(htmlFile, config={'responsive': True})
        htmlFiles[torsionTag] = p.relpath(htmlFile, reporterDir)
    
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


def parse_mol2_file(mol2_file_path):
    """
    Parse a MOL2 file and extract atom and bond information.
    
    Args:
        mol2_file_path (str): Path to the MOL2 file
        
    Returns:
        tuple: (atoms_df, bonds_df) where:
            - atoms_df: DataFrame with columns [atom_id, atom_name, x, y, z, atom_type, charge]
            - bonds_df: DataFrame with columns[bond_id, atom1_id, atom2_id, bond_type]
    """
    with open(mol2_file_path, 'r') as f:
        lines = f.readlines()
    
    atoms_data = []
    bonds_data =[]
    section = None
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('@<TRIPOS>'):
            section = line.replace('@<TRIPOS>', '')
            continue
        
        if not line or section is None:
            continue
        
        if section == 'ATOM':
            parts = line.split()
            if len(parts) >= 9:
                atom_id = int(parts[0])
                atom_name = parts[1]
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                atom_type = parts[5]
                charge = float(parts[8])
                atoms_data.append([atom_id, atom_name, x, y, z, atom_type, charge])
        
        elif section == 'BOND':
            parts = line.split()
            if len(parts) >= 4:
                bond_id = int(parts[0])
                atom1_id = int(parts[1])
                atom2_id = int(parts[2])
                bond_type = parts[3]
                bonds_data.append([bond_id, atom1_id, atom2_id, bond_type])
    
    atoms_df = pd.DataFrame(atoms_data, columns=['atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'charge'])
    bonds_df = pd.DataFrame(bonds_data, columns=['bond_id', 'atom1_id', 'atom2_id', 'bond_type'])
    
    return atoms_df, bonds_df


def build_molecular_graph(atoms_df, bonds_df):
    """
    Build a NetworkX graph from parsed MOL2 data.
    
    Args:
        atoms_df: DataFrame with atom information
        bonds_df: DataFrame with bond information
        
    Returns:
        nx.Graph: Graph with node attributes (atom_name, atom_type, partial_charge, pos) 
                  and edge attributes (bond_length)
    """
    import networkx as nx
    
    G = nx.Graph()
    
    # Create set of valid atom IDs for fast lookup
    valid_atom_ids = set(atoms_df['atom_id'].values)
    
    # Add nodes with attributes
    for _, atom in atoms_df.iterrows():
        G.add_node(
            atom['atom_id'],
            atom_name=atom['atom_name'],
            atom_type=atom['atom_type'],
            partial_charge=atom['charge'],
            pos=(atom['x'], atom['y'], atom['z']),
            x=atom['x'],
            y=atom['y'],
            z=atom['z']
        )
    
    # Add edges with bond lengths (only for bonds where both atoms exist)
    for _, bond in bonds_df.iterrows():
        atom1_id = bond['atom1_id']
        atom2_id = bond['atom2_id']
        
        # Skip bonds where either atom doesn't exist in atoms_df
        if atom1_id not in valid_atom_ids or atom2_id not in valid_atom_ids:
            continue
        
        atom1 = atoms_df[atoms_df['atom_id'] == atom1_id].iloc[0]
        atom2 = atoms_df[atoms_df['atom_id'] == atom2_id].iloc[0]
        
        # Calculate bond length
        dx = atom1['x'] - atom2['x']
        dy = atom1['y'] - atom2['y']
        dz = atom1['z'] - atom2['z']
        bond_length = np.sqrt(dx**2 + dy**2 + dz**2)
        
        G.add_edge(
            atom1_id,
            atom2_id,
            bond_length=bond_length,
            bond_type=bond['bond_type']
        )
    
    return G


def create_interactive_graph_visualization(G, outDir, config):
    """
    Create an interactive 2D molecular graph visualization using plotly.
    Uses custom physics layout for better ChemDraw-style molecular geometry representation.
    
    Args:
        G: NetworkX graph with molecular data
        outDir: Output directory for HTML file
        config: Configuration dictionary
        
    Returns:
        str: Relative path to the generated HTML file
    """
    import plotly.graph_objects as go
    import networkx as nx
    import pandas as pd
    
    reporterDir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    chargeGradientColors = make_charge_gradient_colors()
    
    # Extract data for physics layout mock DataFrames
    atoms_data = [{'atom_id': n, 'x': d['pos'][0], 'y': d['pos'][1], 'z': d['pos'][2]} for n, d in G.nodes(data=True)]
    bonds_data =[{'atom1_id': u, 'atom2_id': v} for u, v in G.edges()]
    
    atoms_df_mock = pd.DataFrame(atoms_data)
    bonds_df_mock = pd.DataFrame(bonds_data)
    
    # Compute 2D layout using the custom physics engine
    pos_2d = compute_physics_layout(atoms_df_mock, bonds_df_mock)
    
    # Extract node attributes
    node_names = []
    node_types = []
    node_charges =[]
    node_elements =[]
    
    for node in G.nodes():
        atom_name = G.nodes[node]['atom_name']
        node_names.append(atom_name)
        node_types.append(G.nodes[node]['atom_type'])
        node_charges.append(G.nodes[node]['partial_charge'])
        
        # Extract element from atom name
        name = atom_name.strip()
        if len(name) >= 2 and name[1].islower():
            element = name[:2]
        else:
            element = name[0]
        node_elements.append(element)
    
    # Define color scheme: C:grey, H:white, O:red, N:blue, S:yellow, Other:orange
    element_colors = {
        'C': '#808080',   # grey
        'H': '#FFFFFF',   # white
        'O': '#FF0000',   # red
        'N': '#0000FF',   # blue
        'S': '#FFFF00',   # yellow
    }
    
    # Create edges using 2D positions
    edge_x = []
    edge_y =[]
    
    for edge in G.edges():
        x0, y0 = pos_2d[edge[0]]
        x1, y1 = pos_2d[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    # Create edge trace
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        mode='lines',
        line=dict(color='#666666', width=3),
        hoverinfo='none',
        showlegend=False
    )
    
    # Map charges to colors using gradient (for background highlight)
    min_charge = min(node_charges)
    max_charge = max(node_charges)
    
    # Normalize charges to [0, 1]
    if max_charge != min_charge:
        normalized_charges =[(c - min_charge) / (max_charge - min_charge) for c in node_charges]
    else:
        normalized_charges = [0.5] * len(node_charges)
    
    # Map normalized charges to gradient colors
    def get_color_from_gradient(value):
        """Map a value in [0,1] to the gradient colors"""
        idx = int(value * (len(chargeGradientColors) - 1))
        idx = max(0, min(idx, len(chargeGradientColors) - 1))
        return chargeGradientColors[idx]
    
    charge_colors = [get_color_from_gradient(nc) for nc in normalized_charges]
    
    # Prepare node data
    node_x = []
    node_y =[]
    node_colors = []
    node_sizes = []
    hover_text =[]
    
    highlight_x = []
    highlight_y = []
    highlight_colors = []
    highlight_sizes =[]
    
    for i, node in enumerate(G.nodes()):
        x, y = pos_2d[node]
        node_x.append(x)
        node_y.append(y)
        
        # Get element for this atom
        element = node_elements[i]
        
        # Determine color based on element
        color = element_colors.get(element, '#FFA500')  # orange for other
        node_colors.append(color)
        
        # Size: hydrogens smaller than heavy atoms
        if element == 'H':
            base_size = 8
            highlight_size = 18
        else:
            base_size = 14
            highlight_size = 28
        
        node_sizes.append(base_size)
        
        # Create highlight circle (larger circle behind with charge color)
        highlight_x.append(x)
        highlight_y.append(y)
        highlight_colors.append(charge_colors[i])
        highlight_sizes.append(highlight_size)
        
        # Hover text
        text = f"Atom: {node_names[i]}<br>"
        text += f"Element: {element}<br>"
        text += f"Type: {node_types[i]}<br>"
        text += f"Charge: {node_charges[i]:.3f}<br>"
        text += f"ID: {node}"
        hover_text.append(text)
    
    # Create highlight trace (larger circles behind atoms showing charge)
    highlight_trace = go.Scatter(
        x=highlight_x,
        y=highlight_y,
        mode='markers',
        marker=dict(
            size=highlight_sizes,
            color=highlight_colors,
            opacity=0.5,
            line=dict(color='#FFFFFF', width=1)
        ),
        hoverinfo='skip',
        showlegend=False
    )
    
    # Create node trace (atoms with element colors)
    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode='markers',
        marker=dict(
            size=node_sizes,
            color=node_colors,
            line=dict(color='#000000', width=1.5)
        ),
        hovertext=hover_text,
        hoverinfo='text',
        showlegend=False
    )
    
    # Create figure with layers: edges, charge highlights, atoms
    fig = go.Figure(data=[edge_trace, highlight_trace, node_trace])
    
    fig.update_layout(
        title=dict(
            text='Molecular Graph Visualization',
            font=dict(color='white', size=20)
        ),
        xaxis=dict(
            visible=False
        ),
        yaxis=dict(
            visible=False,
            scaleanchor='x',
            scaleratio=1
        ),
        paper_bgcolor='black',
        plot_bgcolor='black',
        showlegend=False,
        margin=dict(l=0, r=0, t=40, b=0),
        hovermode='closest'
    )
    
    # Save as HTML
    graph_html_path = p.join(outDir, "molecular_graph.html")
    fig.write_html(graph_html_path)
    
    relative_path = p.relpath(graph_html_path, reporterDir)
    return relative_path
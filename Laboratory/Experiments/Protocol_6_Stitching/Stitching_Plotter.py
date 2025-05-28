import os
from os import path as p
import matplotlib
matplotlib.use('Agg')  # Set the Agg backend before importing pyplot
import matplotlib.pyplot as plt
import numpy as np
from subprocess import call, PIPE
import imageio
from typing import List, Dict # Added for type hinting

# Placeholder classes (extend if needed)
class FilePath:
    pass
class DirectoryPath:
    pass

#####################################################################
def extract_number(filename: str) -> int:
    """Extracts the numerical suffix from a filename (e.g., 'fitting_shuffle_1.png' -> 1)."""
    # Extract the number from the filename
    return int(filename.split('_')[-1].split('.')[0])
    
#####################################################################

def make_gif(in_dir: DirectoryPath, out_gif: FilePath) -> None:
    """
    Creates a GIF from all PNG files in a directory that contain 'fitting_shuffle' in their name.

    Args:
        in_dir: Directory containing the PNG files.
        out_gif: Path to save the output GIF file.
    """
    pngFiles = sorted(
    [file for file in os.listdir(in_dir)  # type: ignore
     if "fitting_shuffle" in file and file.endswith(".png")],
    key=extract_number)    

    images = []
    for file in pngFiles:
        images.append(imageio.v3.imread(p.join(in_dir, file))) # type: ignore

    imageio.mimwrite(out_gif, images, 'GIF', fps=2, loop = 0) # type: ignore



#####################################################################
def set_rc_params() -> None:
    """Sets matplotlib rcParams for consistent plot styling."""
    darkGrey: str = '#1a1a1a'

    yellow: str = '#FFFF00'

    ## set text and axes colors
    plt.rcParams['text.color'] = yellow
    plt.rcParams['axes.labelcolor'] = yellow
    plt.rcParams['xtick.color'] = yellow
    plt.rcParams['ytick.color'] = yellow
    plt.rcParams['axes.edgecolor'] = yellow
    ## set the background color
    plt.rcParams['axes.facecolor'] = darkGrey
    plt.rcParams['figure.facecolor'] = darkGrey
    ## sort out font sizes
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.titlesize'] = 22
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['legend.fontsize'] = 18



#####################################################################
def plot_qmmm_energies(qm_total_energy: np.ndarray,
                        qm_torsion_energy: np.ndarray,
                        mm_total_energy: np.ndarray, 
                        mm_torsion_energy: np.ndarray,
                        cosine_components: Dict[float, np.ndarray],
                        out_dir: DirectoryPath,
                        shuffle_index: int) -> None:
    """
    Plots QM vs MM total energies and QM vs MM torsion energies, including cosine components.

    Args:
        qm_total_energy: NumPy array of QM total energies.
        qm_torsion_energy: NumPy array of QM torsion energies.
        mm_total_energy: NumPy array of MM total energies.
        mm_torsion_energy: NumPy array of MM torsion energies.
        cosine_components: Dictionary mapping frequency to MM cosine component arrays.
        out_dir: Directory to save the plot.
        shuffle_index: Index of the current shuffle/iteration for titling the plot.
    """
    

    set_rc_params()
    ## init some colors to be used
    white :str = '#FFFFFF'
    magenta : str = '#FF00FF'
    brightCyan: str = '#00FFFF'
    ## construct angles
    angles = np.linspace(0, 360, len(qm_total_energy))

    fig, axis = plt.subplots( nrows=1, ncols=2, figsize=(12, 8))  # create figure
    fig.subplots_adjust(bottom=0.25)  # Adjust the bottom margin

    ## plot QM and MM total energies GLOW
    for n in range(1, 7):
        axis[0].plot(angles, qm_total_energy, color=magenta, linewidth=1+n, alpha=0.1)
        axis[0].plot(angles, mm_total_energy, color=brightCyan, linewidth=1+n, alpha=0.1)
    ## plot main traces
    axis[0].plot(angles, qm_total_energy, label='QM', linewidth=2, color=magenta)
    axis[0].plot(angles, mm_total_energy, label='MM', linewidth=2, color=brightCyan)
    ## LEGEND ##
    axis[0].legend(['QM', 'MM'],
                    loc="upper center",
                    bbox_to_anchor=(0.5, -0.12),
                    ncol=2,
                    handlelength=0,
                    handletextpad=0,
                    labelcolor=[magenta, brightCyan])
    
    axis[0].set_xlabel('Torsion Angle')
    axis[0].set_ylabel('Energy (Kcal / mol)')
    axis[0].set_title("Total Energies")
    axis[0].set_xlim(0, 360)
    axis[0].set_xticks(np.arange(0, 361, 60))

    ## plot QM and MM torsion energies GLOW
    for n in range(1, 7):
        axis[1].plot(angles, qm_torsion_energy, color=magenta, linewidth=1+n, alpha=0.1)
        axis[1].plot(angles, mm_torsion_energy, color=brightCyan, linewidth=1+n, alpha=0.1)
    ## plot main traces
    axis[1].plot(angles, qm_torsion_energy, label='QM', linewidth=2, color=magenta)
    axis[1].plot(angles, mm_torsion_energy, label='MM', linewidth=2, color=brightCyan)

    ## plot cosine components
    for freq, cosineComponent in cosine_components.items(): # Renamed variable
        axis[1].plot(angles, cosineComponent, linestyle='--',
                     label=f'Cosine Component {freq:.0f}', color=white, alpha=0.5)
    # Dummy plot for legend
    axis[1].plot([], [], linestyle='--', color=white, alpha=0.5,
                 label='Parameters')
    ## LEGEND ##
    axis[1].legend(['QM', 'MM', 'Parameters'],
                    loc="upper center",
                    bbox_to_anchor=(0.5, -0.12),
                    ncol=3,
                    handlelength=0,
                    handletextpad=0,
                    labelcolor=[magenta, brightCyan, white])

    axis[1].set_xlabel('Torison Angle')
    axis[1].set_ylabel('Energy (Kcal / mol)')
    axis[1].set_title("Torsion Energies")
    axis[1].set_xlim(0, 360)
    axis[1].set_xticks(np.arange(0, 361, 60))

    ## add a shuffle index label
    fig.text(0.1, 0.9, f'{shuffle_index + 1}', fontsize=16,
             color='yellow', bbox=dict(facecolor='none', 
             edgecolor='yellow', boxstyle='round,pad=0.5'))


    fig.savefig(p.join(out_dir, f"fitting_shuffle_{shuffle_index+1}.png")) # type: ignore
    plt.close(fig)
#####################################################################


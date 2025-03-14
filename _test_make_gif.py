import imageio
import os
from os import path as p

def extract_number(filename):
    # Extract the number from the filename
    return int(filename.split('_')[-1].split('.')[0])

def make_gif(inDir, outGif):

    pngFiles = sorted(
    [file for file in os.listdir(inDir) 
     if "fitting_shuffle" in file and file.endswith(".png")],
    key=extract_number)    

    images = []
    for file in pngFiles:
        images.append(imageio.v3.imread(p.join(inDir, file)))

    imageio.mimwrite(outGif, images, 'GIF', fps=1, loop = 0)

make_gif("/home/esp/scriptDevelopment/drFrankenstein/02_NMH_outputs/05_parameter_fitting/qm-mm_parameter_fitting/N-CA-CB-CG",
         "_test.gif")
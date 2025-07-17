from pdbUtils import pdbUtils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--inPdb', '-i', type=str)
parser.add_argument('--outPdb', '-o', type=str)
parser.add_argument('--resName', '-r', type=str)
args = parser.parse_args()


inpPdb = args.inPdb
outPdb = args.outPdb
resName = args.resName
pdbDf = pdbUtils.pdb2df(inpPdb)
pdbDf["RES_NAME"] = resName 
pdbDf["RES_ID"] = 1
pdbUtils.df2pdb(pdbDf, outPdb)
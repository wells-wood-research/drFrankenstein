## BASIC IMPORTS ##
from os import path as p
import os
import pandas as pd
from subprocess import run, PIPE

from pdbUtils import pdbUtils
## CLEAN CODE ##
class FilePath:
    pass
class DirectoryPath:
    pass
from typing import Iterator





def convert_traj_xyz_to_pdb(trajXyzs: list[FilePath],
                           cappedPdb: FilePath,
                           outDir: DirectoryPath) -> Iterator[FilePath]:
    """
    Effectively converts a list of XYZ files to PDB
    In practice, updates a static PDB file with coords from a list of XYZ files

    Args:
        trajXyzs (list[FilePath]): a list of XYZ files
        cappedPdb (FilePath): a PDB file containing correct atom data
        outDir (DirectoryPath): the output dir for this function
    
    Yields:
        FilePath: path to each generated PDB file
    """
    for trajXyz in trajXyzs:
        fileName = p.splitext(p.basename(trajXyz))[0]
        trajPdb = p.join(outDir, f"{fileName}.pdb")
        update_pdb_coords(cappedPdb, trajXyz, trajPdb)
        yield trajPdb

####################
def update_pdb_coords(inPdb: FilePath, xyzFile: FilePath, outPdb: FilePath) -> None:
    """
    updates PDB file with XYZ coords
    
    Args:
        inPdb (FilePath): input PDB file
        xyzFile (FilePath): input XYZ file
        outPdb (FilePath): output PDB file

    Returns:
        None (outPdb already defined!)
    """

    inDf = pdbUtils.pdb2df(inPdb)
    xyzDf = xyz2df(xyzFile)

    inDf["X"] = xyzDf["x"]
    inDf["Y"] = xyzDf["y"]
    inDf["Z"] = xyzDf["z"]

    pdbUtils.df2pdb(inDf, outPdb)

def pdb2xyz(pdbFile: FilePath, xyzFile: FilePath) -> None:
    """
    Uses OpenBabel to convert a PDB file into an XYZ file

    Args:
        pdbFile (FilePath): path to PDB file
        xyzFile (FilePath): path to XYZ file to be created

    Returns:
        None [xyzFile has already been defined]

    """
    obabelCommand = ["obabel", pdbFile, "-O", xyzFile]
    run(obabelCommand, stdout=PIPE, stderr=PIPE)

def pdb2mol2(pdbFile: FilePath, mol2File: FilePath) -> None:
    """
    Uses OpenBabel to convert a PDB file into a MOL2 file

    Args:
        pdbFile (FilePath): path to PDB file
        mol2File (FilePath): path to MOL2 file to be created

    Returns:
        None [mol2File has already been defined]

    """
    obabelCommand = ["obabel", pdbFile, "-O", mol2File]
    run(obabelCommand, stdout=PIPE, stderr=PIPE)
    

def parse_rtf(filePath: str) -> dict:
    """
    Parse a CHARMM RTF file into a structured dictionary.
    
    Args:
        filePath (str): Path to the CHARMM RTF file.
    
    Returns:
        Dict: Parsed RTF data with sections for mass, residues, and other declarations.
    """
    rtfData = {
        "mass": [],  # List of mass entries
        "declarations": [],  # DECL statements
        "defaults": {},  # DEFA statements
        "residues": {},  # RESI definitions
        "auto": []  # AUTO statements
    }
    
    currentResidue = None
    currentGroup = None
    parsingIc = False
    
    with open(filePath, 'r') as file:
        lines = file.readlines()
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith('*'):  # Skip empty lines and comments
            continue
        
        # Parse MASS entries
        if line.startswith('MASS'):
            parts = line.split()
            if len(parts) >= 4:
                massEntry = {
                    "index": parts[1],
                    "atomType": parts[2],
                    "mass": float(parts[3]),
                    "comment": ' '.join(parts[4:]) if len(parts) > 4 else ''
                }
                rtfData["mass"].append(massEntry)
        
        # Parse DECL statements
        elif line.startswith('DECL'):
            rtfData["declarations"].append(line.split()[1])
        
        # Parse DEFA statements
        elif line.startswith('DEFA'):
            parts = line.split()
            rtfData["defaults"] = {
                "first": parts[1] if len(parts) > 1 else '',
                "last": parts[3] if len(parts) > 3 else ''
            }
        
        # Parse AUTO statements
        elif line.startswith('AUTO'):
            rtfData["auto"].extend(line.split()[1:])
        
        # Parse RESI (residue) definitions
        elif line.startswith('RESI'):
            parts = line.split()
            residueName = parts[1]
            charge = float(parts[2]) if len(parts) > 2 else 0.0
            currentResidue = {
                "name": residueName,
                "charge": charge,
                "groups": [],
                "atoms": [],
                "bonds": [],
                "doubleBonds": [],
                "impropers": [],
                "cmap": [],
                "donors": [],
                "acceptors": [],
                "ic": []
            }
            rtfData["residues"][residueName] = currentResidue
        
        # Parse GROUP within residue
        elif line.startswith('GROUP') and currentResidue:
            currentGroup = []
            currentResidue["groups"].append(currentGroup)
        
        # Parse ATOM within residue
        elif line.startswith('ATOM') and currentResidue:
            parts = line.split()
            if len(parts) >= 4:
                atom = {
                    "name": parts[1],
                    "type": parts[2],
                    "charge": float(parts[3]),
                    "comment": ' '.join(parts[4:]) if len(parts) > 4 else ''
                }
                currentResidue["atoms"].append(atom)
                if currentGroup is not None:
                    currentGroup.append(atom["name"])
        
        # Parse BOND
        elif line.startswith('BOND') and currentResidue:
            parts = line.split()[1:]
            for i in range(0, len(parts), 2):
                if i + 1 < len(parts):
                    currentResidue["bonds"].append((parts[i], parts[i + 1]))
        
        # Parse DOUBLE (double bonds)
        elif line.startswith('DOUBLE') and currentResidue:
            parts = line.split()[1:]
            for i in range(0, len(parts), 2):
                if i + 1 < len(parts):
                    currentResidue["doubleBonds"].append((parts[i], parts[i + 1]))
        
        # Parse IMPR (improper dihedrals)
        elif line.startswith('IMPR') and currentResidue:
            parts = line.split()[1:]
            for i in range(0, len(parts), 4):
                if i + 3 < len(parts):
                    currentResidue["impropers"].append(parts[i:i + 4])
        
        # Parse CMAP
        elif line.startswith('CMAP') and currentResidue:
            parts = line.split()[1:]
            currentResidue["cmap"].append(parts)
        
        # Parse DONOR
        elif line.startswith('DONOR') and currentResidue:
            parts = line.split()[1:]
            currentResidue["donors"].append(parts)
        
        # Parse ACCEPTOR
        elif line.startswith('ACCEPTOR') and currentResidue:
            parts = line.split()[1:]
            currentResidue["acceptors"].append(parts)
        
        # Parse IC (internal coordinates)
        elif line.startswith('IC') and currentResidue:
            parts = line.split()
            if len(parts) == 9:
                icEntry = {
                    "atoms": parts[1:5],
                    "star": parts[2] if parts[2].startswith('*') else None,
                    "dist": float(parts[5]),
                    "angle1": float(parts[6]),
                    "dihedral": float(parts[7]),
                    "angle2": float(parts[8]),
                    "dist2": float(parts[9])
                }
                currentResidue["ic"].append(icEntry)
    
    return rtfData



def write_rtf(data, outputFile):
    with open(outputFile, 'w') as f:
        # Write header
        f.write('* Toppar stream file generated by\n')
        f.write('* CHARMM General Force Field (CGenFF) program version 4.0\n')
        f.write('* For use with CGenFF version 4.6\n')
        f.write('*\n\n')
        f.write('read rtf card append\n')
        f.write('* Topologies generated by\n')
        f.write('* CHARMM General Force Field (CGenFF) program version 4.0\n')
        f.write('*\n')
        f.write('36 1\n\n')
        
        # Write penalty comment
        f.write('! "penalty" is the highest penalty score of the associated parameters.\n')
        f.write('! Penalties lower than 10 indicate the analogy is fair; penalties between 10\n')
        f.write('! and 50 mean some basic validation is recommended; penalties higher than\n')
        f.write('! 50 indicate poor analogy and mandate extensive validation/optimization.\n\n')
        
        # Write residue information
        resiName = list(data['residues'].keys())[0]
        resiData = data['residues'][resiName]
        charge = resiData['charge']
        chargePenalty = max(atom['comment'].split()[1] for atom in resiData['atoms'])
        f.write(f'RESI {resiName:<15} {charge:.3f} ! param penalty=   5.000 ; charge penalty= {chargePenalty}\n')
        
        # Write ring comment
        groupAtoms = resiData['groups'][0]
        ringAtoms = [atom for atom in groupAtoms if atom.startswith(('CB', 'CG1', 'CG2', 'CD1', 'CD2', 'CE'))]
        if ringAtoms:
            f.write(f'!RING * 6 {" ".join(ringAtoms)}\n')
        
        # Write GROUP header
        f.write('GROUP            ! CHARGE   CH_PENALTY Z   EL   NB NBE RNG1 TYP RNG2 TYP RNG3 TYP\n')
        
        # Write atoms
        for atom in resiData['atoms']:
            name = atom['name']
            atomType = atom['type']
            charge = atom['charge']
            comment = atom['comment'].split()
            penalty = comment[1]
            z = comment[2]
            el = comment[3]
            nb = comment[4]
            nbe = comment[5]
            extra = ' '.join(comment[6:]) if len(comment) > 6 else ''
            f.write(f'ATOM {name:<7} {atomType:<6} {charge:>6.3f} ! {penalty:>8} {z:>3} {el:<4} {nb:>2} {nbe:<3} {extra}\n')
        
        # Write bonds
        f.write('               ! TYP INR\n')
        for bond in resiData['bonds']:
            if bond[0] == '!':
                continue
            atom1, atom2 = bond
            # Determine bond type (1 for single, 2 for double)
            bondType = '2' if (atom1, atom2) in resiData.get('doubleBonds', []) or (atom2, atom1) in resiData.get('doubleBonds', []) else '1'
            f.write(f'BOND {atom1:<4} {atom2:<4} ! {bondType:<2} 0\n')
        
        # Write impropers
        for improper in resiData['impropers']:
            f.write(f'IMPR {" ".join(f"{atom:<6}" for atom in improper)}\n')
        
        # Write footer
        f.write('\nEND\n')

# Example usage:
# data = your_dictionary
# write_cgenff_topology(data, 'output.top')

def xyz2df(xyzFile: FilePath) -> pd.DataFrame:
    """
    converts an xyz file to a pd.DataFrame

    Args:
        xyzFile (FilePath): input XYZ file

    Returns:
        xyzDf (pd.DataFrame): dataframe containing index, element and coords
    """

    xyzData = []
    atomIndex = 0
    with open(xyzFile, "r") as f:
        lines = f.readlines()[2:]
        for line in lines:
            atomIndex +=1
            lineData = line.split()
            xyzData.append({"index": atomIndex,
                            "element": lineData[0],
                              "x": float(lineData[1]), 
                              "y": float(lineData[2]),
                                "z": float(lineData[3])})
    return pd.DataFrame(xyzData)



def parse_mol2(mol2File):
    atomData = []
    bondData = []
    readingAtoms = False
    readingBonds = False

    with open(mol2File, "r") as mol2:
        for line in mol2:
            if line.strip() == "":
                continue
            if line.startswith("@<TRIPOS>ATOM"):
                readingAtoms=True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                readingBonds=True
                readingAtoms=False
                continue
            if line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                break
            if readingAtoms:
                atomData.append(line.split())
            elif readingBonds:
                bondData.append(line.split())


        
    atomDataColumns = ["ATOM_ID", "ATOM_NAME", "X", "Y", "Z", "ATOM_TYPE", "RES_ID", "RES_NAME", "CHARGE"]
    atomDf = pd.DataFrame(atomData, columns=atomDataColumns)
    bondDataColumns = ["BOND_ID", "ATOM_A_ID", "ATOM_B_ID", "BOND_ORDER"]
    bondDf = pd.DataFrame(bondData, columns=bondDataColumns)

    return atomDf, bondDf
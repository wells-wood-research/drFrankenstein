import os
from os import path as p

import pandas as pd
import numpy as np

from pdbUtils import pdbUtils
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
from itertools import combinations
import math
## CLEAN CODE CLASSES ##
class FilePath:
    pass
class DirectoryPath:
    pass



def dummy_inputs():
    strFile = "/home/esp/scriptDevelopment/drFrankenstein/Inputs/PHG_capped.str"
    cappedPdb = "/home/esp/scriptDevelopment/drFrankenstein/06_PHG_outputs/01_termini_capping/geometry_optimisation/PHG_capped_opt.pdb"
    moleculeName = "PHG"
    outDir = "/home/esp/scriptDevelopment/drFrankenstein/_hidden_scripts/BB_types"
    cgenffRtf = "/home/esp/scriptDevelopment/drFrankenstein/Laboratory/Experiments/Protocol_5_Stitching/CHARMM_protocols/top_all36_cgenff.rtf"
    cgenffPrm = "/home/esp/scriptDevelopment/drFrankenstein/Laboratory/Experiments/Protocol_5_Stitching/CHARMM_protocols/par_all36_cgenff.prm"
    charmmPrm = "/home/esp/scriptDevelopment/drFrankenstein/_hidden_scripts/BB_types/par_all36m_prot.prm"
    
    chargeGroups = {
        "cTerminal": {
            "atoms": ["C", "O"],
            "charge": 0
        },
        "nTerminal": {
            "atoms": ["CA", "N"],
            "charge": 0
        },
        "CB": {
            "atoms": ["CB"],
            "charge": 0
        },
        "ring": {
            "atoms": ["CG1", "CG2", "CD1", "CD2", "CE"],
            "charge": 0
        }
    }

    return strFile, cappedPdb, moleculeName, outDir, cgenffRtf, cgenffPrm, chargeGroups, charmmPrm



def main():
    ## get info that should come from the config
    strFile, cappedPdb, moleculeName, outDir, cgenffRtf, cgenffPrm, chargeGroups, charmmPrm = dummy_inputs()
    ## create dir, already handled!
    os.makedirs(outDir, exist_ok=True)
    ## split stream (.STR) file into RTF and PRM
    moleculeRtf, moleculePrm = split_charmm_str(strFile, moleculeName, outDir)
    ## load cappedPdb into a dataframe
    cappedDf = pdbUtils.pdb2df(cappedPdb)

    ## read RTF into a dict
    parsedRtf = parse_rtf(moleculeRtf)

    ## replace backbone atom atomTypes with backbone defaults
    parsedRtf, cgenffToBackboneDefaultsMap = map_backbone_atoms_back_to_defaults(parsedRtf, moleculeName)


    ## create a PRM dict containing all bonds, angles, dihedrals, impropers
    parsedFullPrm = construct_full_prm(moleculePrm, cgenffPrm, charmmPrm, parsedRtf,  cappedDf, moleculeName)

    parsedCgenffRtf = parse_rtf(cgenffRtf)

    parsedRtf = add_mass_section(parsedRtf, moleculeName, parsedCgenffRtf)


    fullPrm = map_prm_entries_to_backbone_defaults(fullPrm, cgenffToBackboneDefaultsMap)

    parsedRtf = split_atom_groups(parsedRtf, moleculeName, cappedDf, chargeGroups)


    outRtf = p.join(outDir, f"{moleculeName}_edited.rtf")
    write_charmm_rtf(parsedRtf, outRtf)
    write_charmm_prm(fullPrm, p.join(outDir, f"{moleculeName}_edited.prm"))
def map_prm_entries_to_backbone_defaults(parsedPrm, atomTypeMap):

    for interaction, parameters in parsedPrm.items():
        if interaction == "NONBONDED":
            for param in parameters:
                param["atom"] = atomTypeMap.get(param["atom"],param["atom"])
        else:
            for param in parameters:
                param["atoms"] = [atomTypeMap.get(atom, atom) for atom in param["atoms"]]
    return parsedPrm


def construct_full_prm(moleculePrm, cgenffPrm, charmmPrm, parsedRtf,  cappedDf, moleculeName):
    parsedMolPrm = parse_prm(moleculePrm)
    parsedCgenffPrm = parse_prm(cgenffPrm)
    parsedCharmmPrm = parse_prm(charmmPrm)


    atomTypeMap = {atom["name"]: atom["type"] for atom in parsedRtf["residues"][moleculeName]["atoms"]}

    bonds, angles, dihedrals, impropers = find_bond_angle_dihedrals(cappedDf)

    mergedPrm = merge_prms(parsedMolPrm, parsedCgenffPrm, parsedCharmmPrm)

    fullPrm = {}

    fullPrm["BONDS"] = construct_bonds(bonds, atomTypeMap, mergedPrm)
    fullPrm["ANGLES"] = construct_angles(angles, atomTypeMap, mergedPrm)
    fullPrm["DIHEDRALS"] = construct_dihedrals(dihedrals, atomTypeMap, mergedPrm)
    fullPrm["IMPROPERS"] = construct_impropers(impropers, atomTypeMap, mergedPrm)
    fullPrm["NONBONDED"] = construct_nonbonded(parsedCgenffPrm, atomTypeMap)


    return fullPrm

def merge_prms(parsedMolPrm, parsedCgenffPrm, parsedCharmmPrm):

    mergedPrm = {}
    mergedPrm["BONDS"] = parsedMolPrm["BONDS"] + parsedCgenffPrm["BONDS"] + parsedCharmmPrm["BONDS"]
    mergedPrm["ANGLES"] =  parsedMolPrm["ANGLES"] + parsedCgenffPrm["ANGLES"] + parsedCharmmPrm["ANGLES"]
    mergedPrm["DIHEDRALS"] =  parsedMolPrm["DIHEDRALS"] + parsedCgenffPrm["DIHEDRALS"] + parsedCharmmPrm["DIHEDRALS"]
    mergedPrm["IMPROPERS"] =  parsedMolPrm["IMPROPERS"] + parsedCgenffPrm["IMPROPERS"] + parsedCharmmPrm["IMPROPERS"]
    mergedPrm["NONBONDED"] =  parsedMolPrm["NONBONDED"] + parsedCgenffPrm["NONBONDED"] + parsedCharmmPrm["NONBONDED"]

    return mergedPrm

def construct_bonds(bonds, atomTypeMap, mergedPrm):
    """
    Construct a list of unique bond parameters by matching atom types.
    
    Parameters:
    - bonds: List of bond definitions (each a list/tuple of two atom identifiers).
    - atomTypeMap: Dict mapping atom identifiers to their atom types.
    - mergedPrm: Dict containing 'BONDS' key with list of parameter dicts.
                 Each param dict has an 'atoms' key with a list of two atom types.
    
    Returns:
    - List of unique bond parameter dictionaries.
    
    Raises:
    - KeyError: If required keys are missing in atomTypeMap or mergedPrm.
    - ValueError: If a bond has an invalid number of atoms.
    """
    if "BONDS" not in mergedPrm:
        raise KeyError("mergedPrm missing 'BONDS' key")
    
    bondTerms = []
    seenParamIndexes = set()  # Track unique parameters by atom type tuple
    
    for bond in bonds:
        if len(bond) != 2:
            raise ValueError(f"Invalid bond {bond}: must have exactly 2 atoms")
        
        # Map atoms to types
        try:
            atomTypes = [atomTypeMap[atom] for atom in bond]
        except KeyError as e:
            raise KeyError(f"Atom {e} not found in atomTypeMap")
        
        atomTypesReversed = atomTypes[::-1]
        paramsFound = False
        
        for paramIndex, bondParam in enumerate(mergedPrm["BONDS"]):
            paramAtoms = bondParam["atoms"]
            # Check for exact match (forward or reversed)
            if paramAtoms == atomTypes or paramAtoms == atomTypesReversed:
                if paramIndex not in seenParamIndexes:
                    bondTerms.append(bondParam)
                    seenParamIndexes.add(paramIndex)
                paramsFound = True
        
        if not paramsFound:
            print(f"No parameters found for bond with atom types {atomTypes}")

    return bondTerms


def construct_angles(angles, atomTypeMap, mergedPrm):
    """
    Construct a list of unique angle parameters by matching atom types.
    
    Parameters:
    - angles: List of angle definitions (each a list/tuple of three atom identifiers).
    - atomTypeMap: Dict mapping atom identifiers to their atom types.
    - mergedPrm: Dict containing 'ANGLES' key with list of parameter dicts.
                 Each param dict has an 'atoms' key with a list of three atom types.
    
    Returns:
    - List of unique angle parameter dictionaries.
    
    Raises:
    - KeyError: If required keys are missing in atomTypeMap or mergedPrm.
    - ValueError: If an angle has an invalid number of atoms.
    """
    if "ANGLES" not in mergedPrm:
        raise KeyError("mergedPrm missing 'ANGLES' key")
    
    angleTerms = []
    seenParamIndexes = set()  # Track unique parameters by atom type tuple
    
    for angle in angles:
        if len(angle) != 3:
            raise ValueError(f"Invalid angle {angle}: must have exactly 3 atoms")
        
        # Map atoms to types
        try:
            atomTypes = [atomTypeMap[atom] for atom in angle]
        except KeyError as e:
            raise KeyError(f"Atom {e} not found in atomTypeMap")
        
        angleCenterAtom = atomTypes[1]
        angleRadialAtoms = [atomTypes[0], atomTypes[2]]
        angleRadialAtomsReversed = angleRadialAtoms[::-1]
        
        paramsFound = False
        
        for paramIndex, angleParam in enumerate(mergedPrm["ANGLES"]):
            paramAtoms = angleParam["atoms"]
            paramCenterAtom = paramAtoms[1]
            paramRadialAtoms = [paramAtoms[0], paramAtoms[2]]
            
            # Check for exact match (forward or reversed radial atoms)
            if angleCenterAtom == paramCenterAtom:
                if paramRadialAtoms == angleRadialAtoms or paramRadialAtoms == angleRadialAtomsReversed:
                    if paramIndex not in seenParamIndexes:
                        angleTerms.append(angleParam)
                        seenParamIndexes.add(paramIndex)
                    paramsFound = True
        
        if not paramsFound:
            print(f"No parameters found for angle with atom types {atomTypes}")
    
    return angleTerms


def construct_dihedrals(dihedrals, atomTypeMap, mergedPrm):
    """
    Construct a list of unique dihedral parameters by matching atom types.
    
    Parameters:
    - dihedrals: List of dihedral definitions (each a list/tuple of four atom identifiers).
    - atomTypeMap: Dict mapping atom identifiers to their atom types.
    - mergedPrm: Dict containing 'DIHEDRALSphysics: Dict containing 'DIHEDRALS' key with list of parameter dicts.
                 Each param dict has an 'atoms' key with a list of four atom types.
    
    Returns:
    - List of unique dihedral parameter dictionaries.
    
    Raises:
    - KeyError: If required keys are missing in atomTypeMap or mergedPrm.
    - ValueError: If a dihedral has an invalid number of atoms.
    """
    if "DIHEDRALS" not in mergedPrm:
        raise KeyError("mergedPrm missing 'DIHEDRALS' key")
    
    dihedralTerms = []
    seenParamIndexes = set()  # Track unique parameters by atom type tuple
    
    for dihedral in dihedrals:
        if len(dihedral) != 4:
            raise ValueError(f"Invalid dihedral {dihedral}: must have exactly 4 atoms")
        
        # Map atoms to types
        try:
            atomTypes = [atomTypeMap[atom] for atom in dihedral]
        except KeyError as e:
            raise KeyError(f"Atom {e} not found in atomTypeMap")
        
        atomTypesReversed = atomTypes[::-1]
        paramsFound = False
        
        for paramIndex, dihedralParam in enumerate(mergedPrm["DIHEDRALS"]):
            paramAtoms = dihedralParam["atoms"]
            # Check for exact match (forward or reversed)
            if paramAtoms == atomTypes or paramAtoms == atomTypesReversed:
                if paramIndex not in seenParamIndexes:
                    dihedralTerms.append(dihedralParam)
                    seenParamIndexes.add(paramIndex)
                paramsFound = True
        
        if not paramsFound:
            print(f"No parameters found for dihedral with atom types {atomTypes}")
    
    return dihedralTerms



def construct_impropers(impropers, atomTypeMap, mergedPrm):
    """
    Construct a list of unique improper dihedral parameters by matching atom types.
    
    Parameters:
    - impropers: List of improper dihedral definitions (each a list/tuple of four atom identifiers).
    - atomTypeMap: Dict mapping atom identifiers to their atom types.
    - mergedPrm: Dict containing 'IMPROPERS' key with list of parameter dicts.
                 Each param dict has an 'atoms' key with a list of four atom types.
    
    Returns:
    - List of unique improper parameter dictionaries.
    
    Raises:
    - KeyError: If required keys are missing in atomTypeMap or mergedPrm.
    - ValueError: If an improper has an invalid number of atoms.
    """
    if "IMPROPERS" not in mergedPrm:
        raise KeyError("mergedPrm missing 'IMPROPERS' key")
    
    improperTerms = []
    seenParamIndexes = set()  # Track unique parameters by atom type tuple
    
    for improper in impropers:
        if len(improper) != 4:
            raise ValueError(f"Invalid improper {improper}: must have exactly 4 atoms")
        
        # Map atoms to types
        try:
            atomTypes = [atomTypeMap[atom] for atom in improper]
        except KeyError as e:
            raise KeyError(f"Atom {e} not found in atomTypeMap")
        
        improperCenterAtom = atomTypes[0]
        improperRadialAtoms = sorted(atomTypes[1:4]) # Sort for permutation invariance
        paramsFound = False
        
        for paramIndex, improperParam in enumerate(mergedPrm["IMPROPERS"]):
            paramAtoms = improperParam["atoms"]
            paramCenterAtom = paramAtoms[0]
            paramRadialAtoms = sorted(paramAtoms[1:4])
            
            # Check for match (center atom and radial atoms as a set)
            if improperCenterAtom == paramCenterAtom and improperRadialAtoms == paramRadialAtoms:
                if paramIndex not in seenParamIndexes:
                    improperTerms.append(improperParam)
                    seenParamIndexes.add(paramIndex)
                paramsFound = True
        
        if not paramsFound:
            print(f"No parameters found for improper with atom types {atomTypes}")
    
    return improperTerms
def construct_nonbonded(parsedCgenffPrm, atomTypeMap):
    """
    Construct a list of unique nonbonded parameters by matching atom types.
    
    Parameters:
    - parsedCgenffPrm: Dict containing 'NONBONDED' key with list of parameter dicts.
                       Each param dict has an 'atom' key with a single atom type.
    - atomTypeMap: Dict mapping atom identifiers to their atom types.
    
    Returns:
    - List of unique nonbonded parameter dictionaries.
    
    Raises:
    - KeyError: If required keys are missing in parsedCgenffPrm.
    """
    if "NONBONDED" not in parsedCgenffPrm:
        raise KeyError("parsedCgenffPrm missing 'NONBONDED' key")
    
    nonbondedTerms = []
    seenParamIndexes = set()  # Track unique parameters by atom type
    atomTypesForNonbonded = set(atomTypeMap.values())  # Unique atom types
    
    for paramIndex, nonBondedParam in enumerate(parsedCgenffPrm["NONBONDED"]):
        if "atom" not in nonBondedParam:
            raise KeyError(f"Nonbonded parameter at index {paramIndex} missing 'atom' key")
        
        atomType = nonBondedParam["atom"]
        if atomType in atomTypesForNonbonded:
            if paramIndex not in seenParamIndexes:
                nonbondedTerms.append(nonBondedParam)
                seenParamIndexes.add(paramIndex)
    
    return nonbondedTerms


def write_charmm_prm(parsedPrm, outFile):

    with open(outFile, "w") as f:
        if 'BONDS' in parsedPrm:
            f.write('BONDS\n')
            for bond in parsedPrm['BONDS']:
                atom1, atom2 = bond['atoms']
                kValue = bond['k']
                r0Value = bond['r0']
                comment = bond.get('comment', '')
                f.write(f'{atom1:<6} {atom2:<6} {kValue:>8} {r0Value:>8} {comment}\n')
            f.write('\n')

        if 'ANGLES' in parsedPrm:
            f.write('ANGLES\n')
            for angle in parsedPrm['ANGLES']:
                atom1, atom2, atom3 = angle['atoms']
                kValue = angle['k']
                theta0 = angle['theta0']
                comment = angle.get('comment', '')
                f.write(f'{atom1:<6} {atom2:<6} {atom3:<6} {kValue:>8} {theta0:>8} {comment}\n')
            f.write('\n')

        if 'DIHEDRALS' in parsedPrm:
            f.write('DIHEDRALS\n')
            for dihedral in parsedPrm['DIHEDRALS']:
                atom1, atom2, atom3, atom4 = dihedral['atoms']
                kValue = dihedral['k']
                period = dihedral['period']
                phase = dihedral['phase']
                comment = dihedral.get('comment', '')
                f.write(f'{atom1:<6} {atom2:<6} {atom3:<6} {atom4:<6} {kValue:>7} {period:>3} {phase:>8} {comment}\n')
            f.write('\n')

        if 'IMPROPERS' in parsedPrm:
            f.write('IMPROPERS\n')
            for improper in parsedPrm['IMPROPERS']:
                atom1, atom2, atom3, atom4 = improper['atoms']
                kValue = improper['k']
                period = improper['period']
                phase = improper['phase']
                comment = improper.get('comment', '')
                f.write(f'{atom1:<6} {atom2:<6} {atom3:<6} {atom4:<6} {kValue:>7} {period:>3} {phase:>8} {comment}\n')
            f.write('\n')

        if 'NONBONDED' in parsedPrm:
            f.write('NONBONDED\n')
            for nonbonded in parsedPrm['NONBONDED']:
                atom = nonbonded['atom']
                charge = nonbonded['charge']
                epsilon = nonbonded['epsilon']
                ljParam = nonbonded['LJRmin/2']
                comment = nonbonded.get('comment', '')
                f.write(f'{atom:<9} {charge:<10} {epsilon:<12} {ljParam:<7} {comment}\n')
            f.write('\n')
        f.write('END\n')


def create_adjacency_matrix(pdbDf, bondTolerance=1.5):
    """
    Create an adjacency matrix for the molecule based on atomic distances and covalent radii.
    
    Parameters:
    - pdbDf: pandas DataFrame with columns ['ATOM_NAME', 'element', 'x', 'y', 'z']
    - bondTolerance: float, multiplier for covalent radii sum (default: 1.6)
    
    Returns:
    - adj_matrix: pandas DataFrame, adjacency matrix (1 for bonded, 0 otherwise)
    """

    n_atoms = len(pdbDf)
    # Initialize adjacency matrix as a DataFrame of zeros
    adj_matrix = pd.DataFrame(0, index=pdbDf.index, columns=pdbDf.index, dtype=int)
    
    
    # Compute pairwise distances
    for i, j in combinations(range(n_atoms), 2):
        atom_i = pdbDf.iloc[i]
        atom_j = pdbDf.iloc[j]
        # Euclidean distance
        dist = math.sqrt(
            (atom_i['X'] - atom_j['X'])**2 +
            (atom_i['Y'] - atom_j['Y'])**2 +
            (atom_i['Z'] - atom_j['Z'])**2
        )

        if dist <= bondTolerance:
            adj_matrix.iloc[i, j] = 1
            adj_matrix.iloc[j, i] = 1  # Symmetric matrix
    
    return adj_matrix



def _plot_adjacency_matrix(adjacency_matrix, ATOM_NAMEs):
    """
    Plot a heatmap of the adjacency matrix with atom names as labels.
    
    Parameters:
    - adjacency_matrix: pandas DataFrame, adjacency matrix
    - ATOM_NAMEs: pandas Series or list, atom names for labeling
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    # Create a figure
    plt.figure(figsize=(10, 8))
    
    # Create heatmap using seaborn
    sns.heatmap(
        adjacency_matrix,
        xticklabels=ATOM_NAMEs,
        yticklabels=ATOM_NAMEs,
        cmap='binary',  # Black and white color scheme (0: white, 1: black)
        cbar=False,     # No colorbar needed for binary matrix
        square=True     # Make cells square for better visualization
    )
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    # Add labels and title
    plt.xlabel('Atoms')
    plt.ylabel('Atoms')
    plt.title('Molecular Adjacency Matrix')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Show plot
    plt.show()



def find_bond_angle_dihedrals(pdbDf):
    """
    Identify atom names participating in bonds, angles, dihedral angles, and improper angles.
    
    Parameters:
    - pdbDf: pandas DataFrame with columns ['ATOM_NAME', 'element', 'x', 'y', 'z']
    
    Returns:
    - bonds: list of tuples, each containing (ATOM_NAME_i, ATOM_NAME_j)
    - angles: list of tuples, each containing (ATOM_NAME_a, ATOM_NAME_b, ATOM_NAME_c)
    - dihedrals: list of tuples, each containing (ATOM_NAME_a, ATOM_NAME_b, ATOM_NAME_c, ATOM_NAME_d)
    - impropers: list of tuples, each containing (ATOM_NAME_a, ATOM_NAME_b, ATOM_NAME_center, ATOM_NAME_c)
    """
    # Create adjacency matrix (assuming create_adjacency_matrix is defined elsewhere)
    adjacency_matrix = create_adjacency_matrix(pdbDf, bondTolerance=1.5)
    
    # Initialize lists
    bonds = []
    angles = []
    dihedrals = []
    impropers = []
    
    # Extract atom names
    ATOM_NAMEs = pdbDf['ATOM_NAME']
    
    # Find bonds
    for i, j in combinations(range(len(pdbDf)), 2):
        if adjacency_matrix.iloc[i, j] == 1:
            bonds.append((ATOM_NAMEs.iloc[i], ATOM_NAMEs.iloc[j]))
    
    # Find angles
    for b in range(len(pdbDf)):
        # Neighbors of atom b (where adj_matrix[b, :] == 1)
        neighbors = adjacency_matrix.index[adjacency_matrix.iloc[b] == 1].tolist()
        # All pairs of neighbors form an angle with b as the central atom
        for a, c in combinations(neighbors, 2):
            angles.append((ATOM_NAMEs.loc[a], ATOM_NAMEs.loc[b], ATOM_NAMEs.loc[c]))
    
    # Find dihedrals
    for b, c in bonds:
        b_idx = pdbDf.index[pdbDf['ATOM_NAME'] == b].tolist()[0]
        c_idx = pdbDf.index[pdbDf['ATOM_NAME'] == c].tolist()[0]
        # Neighbors of b (excluding c)
        neighbors_b = adjacency_matrix.index[
            (adjacency_matrix.loc[b_idx] == 1) & (adjacency_matrix.index != c_idx)
        ].tolist()
        # Neighbors of c (excluding b)
        neighbors_c = adjacency_matrix.index[
            (adjacency_matrix.loc[c_idx] == 1) & (adjacency_matrix.index != b_idx)
        ].tolist()
        for a in neighbors_b:
            for d in neighbors_c:
                dihedrals.append((
                    ATOM_NAMEs.loc[a],
                    ATOM_NAMEs.loc[b_idx],
                    ATOM_NAMEs.loc[c_idx],
                    ATOM_NAMEs.loc[d]
                ))
    
    # Find improper angles
    for center in range(len(pdbDf)):
        # Neighbors of the central atom
        neighbors = adjacency_matrix.index[adjacency_matrix.iloc[center] == 1].tolist()
        # If the central atom has at least 3 neighbors, form improper angles
        if len(neighbors) >= 3:
            # All combinations of 3 neighbors
            for a, b, c in combinations(neighbors, 3):
                impropers.append((
                    ATOM_NAMEs.loc[a],
                    ATOM_NAMEs.loc[b],
                    ATOM_NAMEs.loc[center],  # Central atom as third index
                    ATOM_NAMEs.loc[c]
                ))
    
    return sorted(bonds), sorted(angles), sorted(dihedrals), sorted(impropers)

def split_atom_groups(parsedRtf, moleculeName, cappedDf, chargeGroups):
    newGroups = []
    ## remove capping groups
    uncappedDf = cappedDf[~cappedDf["RES_NAME"].isin(["NME", "ACE"])]

    ## deal with user-defined charge groups
    userDefinedAtoms = []

    for chargeGroupName, chargeGroupData in chargeGroups.items():
        ## get heavy atom names
        heavyChargeGroupAtoms = chargeGroupData["atoms"]
        ## fill in hydrogen atom names
        protonNames = []
        for heavyAtom in heavyChargeGroupAtoms:
            boundProtons = find_bonded_atoms(cappedDf, heavyAtom)
            protonNames.extend(boundProtons)
        protonNames = [atomName for atomName in protonNames if atomName.startswith("H")] ##TODO: make this nicer
        chargeGroupAtoms = heavyChargeGroupAtoms + protonNames
        newGroups.append(chargeGroupAtoms)


    ## deal with left-over atoms, these will all go in one group
    leftOverDf = uncappedDf[~uncappedDf["ATOM_NAME"].isin(userDefinedAtoms)]
    leftOverAtoms = leftOverDf["ATOM_NAME"].to_list()
    if len(leftOverAtoms) > 0:
        newGroups.append(leftOverAtoms)
    
 
    ## deal with Terminal Caps 
    for capResName in ["NME", "ACE"]:
        capDf = cappedDf[cappedDf["RES_NAME"]==capResName]
        for capIndex, capResDf in capDf.groupby("RES_ID"):
            chargeGroupAtoms = capResDf["ATOM_NAME"].to_list() 
            newGroups.append(chargeGroupAtoms)
    parsedRtf["residues"][moleculeName]["groups"] = newGroups
    return parsedRtf
      


def calculate_distance(row: pd.Series,
                        atomX: float,
                          atomY: float,
                            atomZ: float):
    """
    Calculates distance between a row of a pdb dataframe and a set of coords
    Used in an apply function
    
    Args:
        row (pd.Series): row of a pdb dataframe
        atomX (float): x coord
        atomY (float): y coord
        atomZ (float): z coord
    Returns:
        distance (float): distance
    """
    return np.sqrt((row['X'] - atomX) ** 2 +
                   (row['Y'] - atomY) ** 2 +
                   (row['Z'] - atomZ) ** 2)

def show_dict(parsedRtf):
    for key, value in parsedRtf.items():
        if type(value) == dict:
            for key2, value2 in value.items():
                print(key2, value2)
        elif type(value) == list:
            print(key)
            for value2 in value:
                print(value2)
        else:
            print(key,value)



def add_mass_section(parsedRtf, moleculeName, parsedCgenffRtf):
    atomTypes = list(set([atom["type"] for atom in parsedRtf["residues"][moleculeName]["atoms"]]))
    massSection = []

    for massEntry in parsedCgenffRtf["mass"]:
        if massEntry["atomType"] in atomTypes:
            massSection.append(massEntry)
    parsedRtf["mass"] = massSection
    return parsedRtf

def map_backbone_atoms_back_to_defaults(parsedRtf, moleculeName):
    charmmDefaultMap = {
                    "N" : "NH1",
                    "HN": "H",
                    "CA": "C",
                    "HA": "HB1",
                    "C" : "C",
                    "O" : "O",
                    "NN": "NH1",
                    "HNN1": "H",
                    "CN": "C", 
                    "HCN1": "HB1",
                    "HCN2": "HB1",
                    "HCN3": "HB1",
                    "CC1": "C",
                    "OC": "O",
                    "CC2": "CT1", 
                    "HC1": "HB1",
                    "HC2": "HB1",
                    "HC3": "HB1"
    }

    cgennffToCharmmMap = {}

    for atom in parsedRtf["residues"][moleculeName]["atoms"]:
        if atom["name"] in charmmDefaultMap:
            cgennffToCharmmMap[atom["type"]] = charmmDefaultMap[atom["name"]]

            atom["type"] = charmmDefaultMap[atom["name"]]

    for massEntry in parsedRtf["mass"]:
        if massEntry["atomType"] in cgennffToCharmmMap:
            massEntry["atomType"] =  cgennffToCharmmMap[massEntry["atomType"]]
    ## remove duplicated MASS entries
    parsedRtf["mass"] = [d for d in dict((entry["atomType"], entry) for entry in parsedRtf["mass"]).values()]


    return parsedRtf, cgennffToCharmmMap






def split_charmm_str(strFile, moleculeName, outDir) -> dict:
    """
    Splits a CHARMM STR (stream) file into RTF and PRM files

    Args: 
        config (dict): contains all run info
    Returns:
        config (dict): updated config 
    """

    # moleculeName = config["moleculeInfo"]["moleculeName"]
    # cappedStr = config["runtimeInfo"]["madeByCapping"]["cappedMoleculeStr"]
    # mmTorsionCalculationDir = config["runtimeInfo"]["madeByStitching"]["mmTorsionCalculationDir"]


    moleculePrm = p.join(outDir, f"{moleculeName}.prm")
    moleculeRtf = p.join(outDir, f"{moleculeName}.rtf")

    with open (strFile, "r") as stream, open(moleculePrm, "w") as prm, open(moleculeRtf, "w") as rtf:
        writeRtf = True
        writePrm = False
        for line in stream:
            if line.startswith("read param card flex append"):
                writeRtf = False
                writePrm = True
            if writeRtf:
                rtf.write(line)
            if writePrm:
                prm.write(line)


    return moleculeRtf, moleculePrm


import re

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
    
    with open(filePath, 'r') as file:
        lines = file.readlines()
    
    # Regex pattern for extracting valid floating-point numbers (e.g., -12.34, 0.0, 123)
    float_pattern = r'-?\d*\.?\d+'
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith('*'):  # Skip empty lines and comments
            continue
        
        # Parse MASS entries
        if line.startswith('MASS'):
            parts = line.split()
            if len(parts) >= 4:
                # Extract numeric part using regex
                mass_match = re.search(float_pattern, parts[3])
                mass_value = float(mass_match.group()) if mass_match else 0.0
                massEntry = {
                    "index": parts[1],
                    "atomType": parts[2],
                    "mass": mass_value,
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
            # Extract numeric part using regex
            charge_match = re.search(float_pattern, parts[2]) if len(parts) > 2 else None
            charge = float(charge_match.group()) if charge_match else 0.0
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
                # Extract numeric part using regex
                charge_match = re.search(float_pattern, parts[3])
                atom_charge = float(charge_match.group()) if charge_match else 0.0
                atom = {
                    "name": parts[1],
                    "type": parts[2],
                    "charge": atom_charge,
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
                # Extract numeric parts using regex
                dist_match = re.search(float_pattern, parts[5])
                angle1_match = re.search(float_pattern, parts[6])
                dihedral_match = re.search(float_pattern, parts[7])
                angle2_match = re.search(float_pattern, parts[8])
                dist2_match = re.search(float_pattern, parts[9])
                
                icEntry = {
                    "atoms": parts[1:5],
                    "star": parts[2] if parts[2].startswith('*') else None,
                    "dist": float(dist_match.group()) if dist_match else 0.0,
                    "angle1": float(angle1_match.group()) if angle1_match else 0.0,
                    "dihedral": float(dihedral_match.group()) if dihedral_match else 0.0,
                    "angle2": float(angle2_match.group()) if angle2_match else 0.0,
                    "dist2": float(dist2_match.group()) if dist2_match else 0.0
                }
                currentResidue["ic"].append(icEntry)
    
    return rtfData
def find_bonded_atoms(molDf: pd.DataFrame,
                       atomName: str,
                         distanceCutoff: float = 1.5) -> list:
    """
    Finds all atoms within a certain distance of a given atom

    Args:
        molDf (pd.DataFrame): dataframe of the molecule
        atomName (str): name of the atom to find bonded atoms for
        distanceCutoff (float): distance cutoff
    
    Returns:
        bondedAtoms (list): list of bonded atoms
    
    """
    ## make a copy of input dataframe
    tmpDf = molDf.copy()
    ## make a dataframe containing only our atom of interest
    atomDf = molDf[molDf["ATOM_NAME"] == atomName]
    ## get atom coords
    atomX, atomY, atomZ = atomDf.loc[:,['X', 'Y', 'Z']].astype(float).iloc[0]
    ## calculate distance to all other atoms
    tmpDf['distance_to_atom'] = tmpDf.apply(calculate_distance,
                                             axis=1, atomX=atomX, atomY=atomY, atomZ=atomZ)
    ## find all atoms within distance cutoff
    bondedDf = tmpDf[tmpDf['distance_to_atom'] < distanceCutoff]
    ## convert to list
    bondedAtoms = bondedDf["ATOM_NAME"].to_list()
    ## remove atom of interest
    bondedAtoms = [atom for atom in bondedAtoms if atom != atomName]
    return bondedAtoms


def write_charmm_rtf(parsedRtf: dict, outRtf: str) -> None:
    """
    Write a parsed RTF dictionary back to a CHARMM RTF file.
    
    Args:
        parsedRtf (dict): Parsed RTF data dictionary.
        outRtf (str): Path to the output RTF file.
    """
    with open(outRtf, 'w') as file:
        # Write header (optional, can be customized)
        file.write("* CHARMM RTF file edited by drFRANKENSTEIN\n")
        file.write("* \n")
        
        # Write MASS entries
        if parsedRtf.get("mass"):
            for massEntry in parsedRtf["mass"]:
                commentText = f"  {massEntry['comment']}" if massEntry["comment"] else ""
                file.write(
                    f"MASS {massEntry['index']:>4} {massEntry['atomType']:<6} "
                    f"{massEntry['mass']:>8.6f}{commentText}\n"
                )
        
        # Write DECL entries
        if parsedRtf.get("declarations"):
            file.write("\nDECL\n")
            for declEntry in parsedRtf["declarations"]:
                file.write(f"DECL {declEntry}\n")
        
        # Write DEFA entries
        if parsedRtf.get("defaults"):
            defaultSettings = parsedRtf["defaults"]
            if defaultSettings.get("first") or defaultSettings.get("last"):
                file.write(f"\nDEFA {defaultSettings['first']} {defaultSettings['last']}\n")
        
        # Write AUTO entries
        if parsedRtf.get("auto"):
            file.write(f"\nAUTO {' '.join(parsedRtf['auto'])}\n")
        
        # Write RESI entries
        if parsedRtf.get("residues"):
            for resiName, resiData in parsedRtf["residues"].items():
                file.write(f"\nRESI {resiName} {resiData['charge']:.6f}\n")
                
                # Write GROUP entries
                for groupAtoms in resiData.get("groups", []):
                    file.write("GROUP\n")
                    # Write ATOM entries for this group
                    for atomName in groupAtoms:
                        for atomEntry in resiData["atoms"]:
                            if atomEntry["name"] == atomName:
                                commentText = f"  {atomEntry['comment']}" if atomEntry["comment"] else ""
                                file.write(
                                    f"ATOM {atomEntry['name']:<4} {atomEntry['type']:<6} "
                                    f"{atomEntry['charge']:>8.6f}{commentText}\n"
                                )
                
                # Write BOND entries
                if resiData.get("bonds"):
                    for bondPair in resiData["bonds"]:
                        # Skip comment-like entries (e.g., ('!', '1'))
                        if bondPair[0].startswith('!'):
                            continue
                        file.write(f"BOND {bondPair[0]:<4} {bondPair[1]:<4}\n")
                
                # Write DOUBLE entries
                if resiData.get("doubleBonds"):
                    file.write("DOUBLE")
                    for doubleBond in resiData["doubleBonds"]:
                        file.write(f" {doubleBond[0]} {doubleBond[1]}\n")
                    file.write("\n")
                
                # Write IMPR entries
                if resiData.get("impropers"):
                    for improperEntry in resiData["impropers"]:
                        file.write(f"IMPR {improperEntry[0]} {improperEntry[1]} {improperEntry[2]} {improperEntry[3]}\n")
                
                # Write CMAP entries
                if resiData.get("cmap"):
                    file.write("CMAP")
                    for cmapEntry in resiData["cmap"]:
                        file.write(f" {' '.join(cmapEntry)}")
                    file.write("\n")
                
                # Write DONOR entries
                if resiData.get("donors"):
                    file.write("DONOR")
                    for donorEntry in resiData["donors"]:
                        file.write(f" {' '.join(donorEntry)}")
                    file.write("\n")
                
                # Write ACCEPTOR entries
                if resiData.get("acceptors"):
                    file.write("ACCEPTOR")
                    for acceptorEntry in resiData["acceptors"]:
                        file.write(f" {' '.join(acceptorEntry)}")
                    file.write("\n")
                
                # Write IC entries
                if resiData.get("ic"):
                    file.write("IC\n")
                    for icEntry in resiData["ic"]:
                        atomList = icEntry["atoms"]
                        starMark = icEntry["star"] if icEntry["star"] else atomList[1]
                        file.write(
                            f"IC {atomList[0]:<4} {starMark:<4} {atomList[2]:<4} {atomList[3]:<4} "
                            f"{icEntry['dist']:>8.4f} {icEntry['angle1']:>8.4f} {icEntry['dihedral']:>8.4f} "
                            f"{icEntry['angle2']:>8.4f} {icEntry['dist2']:>8.4f}\n"
                        )
        
        # Write end of file
        file.write("\nEND\n")


def parse_prm(moleculePrm: FilePath) -> dict:
    """
    Reads through a CHARMM 
    """
    ## init a bunch of bools
    readingBonds: bool      = False
    readingAngles: bool     = False
    readingDihedrals: bool  = False
    readingImpropers: bool  = False
    readingNonBonded: bool  = False

    ## init empty parsedPrm dict
    parsedPrm = {"BONDS": [], "ANGLES": [], "DIHEDRALS": [], "IMPROPERS": [], "NONBONDED":[]}

    with open(moleculePrm, "r") as f:
        for line in f.readlines():
            if line.strip() == "":
                continue
            elif line.strip().startswith("!"):
                continue
            elif line.startswith("cutnb"):
                continue
            if line.startswith("BONDS"):
                readingBonds = True
            elif line.startswith("ANGLES"):
                readingBonds = False
                readingAngles = True
            elif line.startswith("DIHEDRALS"):
                readingAngles = False
                readingDihedrals = True
            elif line.startswith("IMPROPER"):
                readingDihedrals = False
                readingImpropers = True
            elif line.startswith("NONBONDED"):
                readingNonBonded = True
                readingImpropers = False
            elif line.startswith("END") or line.startswith("NBFIX"):
                break
            else:
                lineData = line.split()

                if "!" in line:
                    comment = "!" + line.split("!")[1].strip()
                else: comment = ""
                if readingBonds:
                    lineParsed = {"atoms": [lineData[0], lineData[1]],
                                    "k": lineData[2], "r0": lineData[3],
                                    "comment": comment}
                    parsedPrm["BONDS"].append(lineParsed)
                elif readingAngles:
                    lineParsed = {"atoms": [lineData[0], lineData[1], lineData[2]],
                                    "k": lineData[3], "theta0": lineData[3],
                                    "comment": comment}
                    parsedPrm["ANGLES"].append(lineParsed)
                elif readingDihedrals:

                    lineParsed = {"atoms": [lineData[0], lineData[1], lineData[2], lineData[3]],
                                    "k": lineData[4], "period": lineData[5], "phase": lineData[6],
                                    "comment": comment}
                    parsedPrm["DIHEDRALS"].append(lineParsed)
                elif readingImpropers:
                    print(line)

                    lineParsed = {"atoms": [lineData[0], lineData[1], lineData[2], lineData[3]],
                                    "k": lineData[4], "period": lineData[5], "phase": lineData[6],
                                    "comment": comment}
                    parsedPrm["IMPROPERS"].append(lineParsed)
                elif readingNonBonded:
                    lineParsed = {"atom" : lineData[0], "charge": lineData[1],
                                   "epsilon": lineData[2], "LJRmin/2": lineData[3],
                                   "comment": comment}
                    parsedPrm["NONBONDED"].append(lineParsed)

    return parsedPrm

if __name__ == "__main__":
    main()
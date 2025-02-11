## BASIC IMPORTS ##
import os
from os import path as p


## RDKIT LIBRARIES ##
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles


## CUSTOM LIBRARIES ##
from pdbUtils import pdbUtils
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def dummy_inputs() -> dict:
    config: dict = {
        "pathInfo": {
            "inputDir": "/home/esp/scriptDevelopment/drFrankenstein/inputs",
            "outputDir": "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs"
        },
        "moleculeInfo": {
            "charge": 0,
            "multiplicity": 1, 
            "moleculePdb": "/home/esp/scriptDevelopment/drFrankenstein/inputs/ALA.pdb"
        }, 
        "scanInfo": {
            "nScanSteps": 37,
            "convergenceTolerance": 0.1,
            "minBatches": 5,
            "maxBatches": 30
        },
        "hardwareInfo": {
            "nCores": 16
        }
    }
    return config
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def capping_protocol(config: dict):

    outputDir = config["pathInfo"]["outputDir"]
    molPdb = config["moleculeInfo"]["moleculePdb"]
    

    inMol = Chem.MolFromPDBFile(molPdb, removeHs=True)
    cappedMol, capIndexNameMap = add_capping_groups(inMol)

    cappedMol = enforce_carboxyl_bonding(cappedMol)

    Chem.SanitizeMol(cappedMol)
    cappedMol = Chem.AddHs(cappedMol)
    AllChem.EmbedMolecule(cappedMol)

    cappingDir = p.join(outputDir, "01_capped_amino_acids")
    os.makedirs(cappingDir,exist_ok=True)

    molName = p.basename(molPdb).split(".")[0]
    cappedMolPdb = p.join(cappingDir, f"{molName}_capped.pdb")
    
    write_pdb(cappedMol, cappedMolPdb)

    rename_capping_atoms(cappedMolPdb, capIndexNameMap)
    return cappedMolPdb
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def rename_capping_atoms(cappedMolPdb, capIndexNameMap):
    cappedDf = pdbUtils.pdb2df(cappedMolPdb)

    cappedDf.loc[cappedDf['ATOM_ID'].isin(capIndexNameMap.keys()),
                  'ATOM_NAME'] = cappedDf['ATOM_ID'].map(capIndexNameMap)
    
    # cappedDf.loc[cappedDf['ATOM_NAME'].isin(["NN", "CN"]), "RES_NAME"] = "NME"
    # cappedDf.loc[cappedDf['ATOM_NAME'].isin(["CC1", "CC2", "OC"]), "RES_NAME"] = "ACE"

    cappedDf["RES_NAME"] = "MOL"
    cappedDf["RES_ID"] = 1

    pdbUtils.df2pdb(cappedDf, cappedMolPdb)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def write_pdb(mol, filename):
    with rdmolfiles.PDBWriter(filename) as writer:
        writer.write(mol)
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def print_smarts(mol):
    smarts = Chem.MolToSmarts(mol)
    print("SMARTS:", smarts)

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def enforce_carboxyl_bonding(mol):
    # Define the initial and target SMARTS patterns
    initial_smarts = '[#7]-[#6@](-[*]([#8])[#7])-[#6]'
    
    # Find the substructure match
    initial_pattern = Chem.MolFromSmarts(initial_smarts)
    matches = mol.GetSubstructMatches(initial_pattern)
    
    if not matches:
        raise ValueError("No matching substructure found.")
    
    # Modify the first match found
    emol = Chem.EditableMol(mol)
    for match in matches:
        # Get the indices of the atoms in the match
        carbon_idx = match[2]
        oxygen_idx = match[3]
        
        # Remove the single bond and add a double bond
        emol.RemoveBond(carbon_idx, oxygen_idx)
        emol.AddBond(carbon_idx, oxygen_idx, Chem.BondType.DOUBLE)
    
    # Return the modified molecule
    return emol.GetMol()
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def add_capping_groups(mol):
    # Use SMARTS to identify the backbone N-terminal nitrogen atom
    backbonePattern = Chem.MolFromSmarts('[#7]-[#6@](-[*]([#8])[#8])-[#6]')
    
    backboneMatches = mol.GetSubstructMatches(backbonePattern)
    if backboneMatches:
        nMethyCappedMol, nMethylIndexNameMap = add_nMethyl_cap(mol, backboneMatches)
        bothCappedMol, acetylIndexNameMap = add_acetyl_cap(nMethyCappedMol, backboneMatches)
    else:
        raise ValueError("No matching substructure found.")
    capIndexNameMap = {**acetylIndexNameMap, **nMethylIndexNameMap}
    return bothCappedMol, capIndexNameMap
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def add_acetyl_cap(mol, backboneMatches):
    n_terminal_idx = backboneMatches[0][0]
    n_terminal = mol.GetAtomWithIdx(n_terminal_idx)
    
    if n_terminal is None:
        raise ValueError("N-terminal nitrogen not found.")
    
    # Add acetyl group (CH3CO)
    acetyl = Chem.MolFromSmiles('C(=O)C')
    mol = Chem.CombineMols(mol, acetyl)
    emol = Chem.EditableMol(mol)
    
    # Connect the acetyl group to the N-terminal
    n_idx = n_terminal.GetIdx()
    c_idx = mol.GetNumAtoms() - 3  # Assuming acetyl is added at the end
    emol.AddBond(n_idx, c_idx, Chem.BondType.SINGLE)
    ## update mol
    mol = emol.GetMol()
    acetylAtomNames = ["CAC1", "CAC2", "OAC"]
    nAtoms = mol.GetNumAtoms()
    acetylIndexNameMap = {nAtoms: "CC2", nAtoms -1 : "OC", nAtoms - 2: "CC1"}

    # Assign names to the acetyl atoms

    return emol.GetMol(), acetylIndexNameMap
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def add_nMethyl_cap(mol, backboneMatches):
    # Identify the C-terminal carbon and oxygen
    c_terminal_carbon_idx = backboneMatches[0][-4]
    c_terminal_oxygen_idx = backboneMatches[0][-2]
    
    # Create an editable molecule
    emol = Chem.EditableMol(mol)
    
    # Remove the C-terminal oxygen
    emol.RemoveAtom(c_terminal_oxygen_idx)
    
    # Add N-methyl group (NCH3)
    nMethyl = Chem.MolFromSmiles('NC')
    mol_with_nMethyl = Chem.CombineMols(emol.GetMol(), nMethyl)
    emol = Chem.EditableMol(mol_with_nMethyl)
    
    # Connect the N-methyl group to the C-terminal carbon
    n_idx = mol_with_nMethyl.GetNumAtoms() - 2  # Assuming N-methyl is added at the end
    emol.AddBond(c_terminal_carbon_idx, n_idx, Chem.BondType.SINGLE)
    ## update mol
    mol = emol.GetMol()
    nAtoms = mol.GetNumAtoms()
    nMethylIndexNameMap = {nAtoms : "CN", nAtoms - 1 : "NN"}

    
    return emol.GetMol(), nMethylIndexNameMap
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
if __name__ == "__main__":
    capping_protocol()
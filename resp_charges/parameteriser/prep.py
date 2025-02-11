import os
from pathlib import Path

# import parmed
from rdkit import Chem
from rdkit.Chem import AllChem

import parameteriser.check_structure as check_structure


def build_molecule_from_smiles(smiles):
    # mol = Chem.MolFromSmiles(smiles)
    mol, res_type, backbone_list, capping_list, sidechain_list = check_structure.check_smiles(smiles)
    #print(capping_list)

    AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=-1) # number of conformers determined by consecutive runs comparison on S-pentenyl and R-pentenyl glycine
    rmslist = []
    AllChem.AlignMolConformers(mol, RMSlist=rmslist)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=1000) # include minimisation iterations to make sure all conformers have been minimised

    return mol, res_type, backbone_list, capping_list, sidechain_list


""""
def build_molecule_from_side_chain_smiles(smiles, stereo):
    backbone = "pass"
"""


def write_mol_to_pdb(mol: Chem, outputfile: Path) -> None:
    Chem.MolToPDBFile(mol, outputfile)


# def get_mol_coordinates(mol):
#     struct = parmed.rdkit.load_rdkit(mol)
#     return struct.positions


def get_dih_random(mol):
    AllChem.SetDihedralDeg(mol.GetConformer(0), 2, 3, 4, 5, randrange(360))
    AllChem.SetDihedralDeg(mol.GetConformer(0), 3, 4, 5, 6, randrange(360))
    AllChem.SetDihedralDeg(mol.GetConformer(0), 0, 1, 2, 3, randrange(360))
    AllChem.SetDihedralDeg(mol.GetConformer(0), 7, 1, 2, 3, randrange(360))

    return mol


def parameterise_mol(mol2_file, ffadd):
    os.system(f"parmchk2 -i  {mol2_file}  -f mol2 -o {ffadd}  -frc  ff14SB -s 2 -a Y")


if __name__ == "__main__":
    smile = "C[C@](CCCC=C)(C(=O)NC)NC(=O)C"
    mol, res_type, backbone_list, capping_list = build_molecule_from_smiles(smile)
    print("capping atoms :" + ",".join(map(str, capping_list)))
    write_mol_to_pdb(mol, outputfile="try.pdb")

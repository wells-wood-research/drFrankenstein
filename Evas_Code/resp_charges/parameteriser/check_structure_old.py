from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from parameteriser.run_resp import run_resp


def rename_atoms(mol):
    return mol


def get_match_bond_indices(mol, match_atom_indices):
    bond_indices = []
    frag_types = []
    for query_bond in mol.GetBonds():
        if (
            query_bond.GetBeginAtomIdx() in match_atom_indices
            and query_bond.GetEndAtomIdx() not in match_atom_indices
        ):
            bond_indices.append(query_bond.GetIdx())
            at = query_bond.GetBeginAtomIdx()
            frag_types.append(get_type(mol, at))

        elif (
            query_bond.GetBeginAtomIdx() not in match_atom_indices
            and query_bond.GetEndAtomIdx() in match_atom_indices
        ):
            bond_indices.append(query_bond.GetIdx())
            at = query_bond.GetEndAtomIdx()
            frag_types.append(get_type(mol, at))
    return bond_indices, frag_types


def get_type(mol, at):
    if (
        mol.GetAtomWithIdx(at).GetSymbol() == "C"
        and str(mol.GetAtomWithIdx(at).GetHybridization()) == "SP2"
    ):
        print(mol.GetAtomWithIdx(at).GetHybridization())
        return "capping"
    elif mol.GetAtomWithIdx(at).GetSymbol() == "N":
        return "capping"
    elif (
        mol.GetAtomWithIdx(at).GetSymbol() == "C"
        and str(mol.GetAtomWithIdx(at).GetHybridization()) == "SP3"
    ):
        return "sidechain"
    else:
        raise Exception("Probleme parsing the structure")


def check_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    patt_backbone = Chem.MolFromSmarts("C(C=O)N([H])")
    mol = Chem.AddHs(mol)

    if mol.HasSubstructMatch(patt_backbone) is False:
        raise Exception("Could not recognise backbone pattern in your molecule")

    match (len(mol.GetSubstructMatches(patt_backbone))):
        case 1:
            print("Treating your molecule as a unstapled residue")

            atom_indices = mol.GetSubstructMatch(patt_backbone)
            bond_indices, frag_types = get_match_bond_indices(mol, atom_indices)
            mol1_f = Chem.FragmentOnBonds(mol, bond_indices, addDummies=False)
            fragments_mols = Chem.GetMolFrags(mol1_f, sanitizeFrags=False)
            frag_types.insert(
                fragments_mols.index(tuple(sorted(atom_indices))), "backbone"
            )
            backbone, capping_groups, sidechain = [], [], []

            for t, frag in zip(frag_types, fragments_mols):
                match t:
                    case "capping":
                        capping_groups.extend(frag)
                    case "backbone":
                        backbone.extend(frag)
                    case "sidechain":
                        sidechain.extend(frag)
            return (mol, 1, backbone, capping_groups)
        case 2:
            print("Treating your residue as a stapled residue")

            atom_indices_1 = mol.GetSubstructMatches(patt_backbone)[0]
            atom_indices_2 = mol.GetSubstructMatches(patt_backbone)[1]

            bond_indices, frag_types = get_match_bond_indices(
                mol, atom_indices_1 + atom_indices_2
            )
            mol1_f = Chem.FragmentOnBonds(mol, bond_indices, addDummies=False)
            fragments_mols = Chem.GetMolFrags(mol1_f, sanitizeFrags=False)

            print(fragments_mols)
            print(sorted(atom_indices_1))
            frag_types.insert(
                fragments_mols.index(tuple(sorted(atom_indices_1))), "backbone"
            )
            frag_types.insert(
                fragments_mols.index(tuple(sorted(atom_indices_2))), "backbone"
            )
            backbone, capping_groups, sidechain = [], [], []

            for t, frag in zip(frag_types, fragments_mols):
                match t:
                    case "capping":
                        capping_groups.extend(frag)
                    case "backbone":
                        backbone.extend(frag)
                    case "sidechain":
                        sidechain.extend(frag)

            return (mol, 2, backbone, capping_groups)

        case _:
            raise ExceptionError(
                " 3 or more backbones fragment found in your molecule. "
                "Your smile chain might be alright but that case in not accepted in current code."
                " please contact developper to include it"
            )
            return 0


if __name__ == "__main__":
    smiles = "C[C@](CCCC=C)(C(=O)NC)NC(=O)C"
    print(check_smiles(smiles))
    smiles = "CNC(CCCC(NC)C(C)=O)C(C)=O"
    print(check_smiles(smiles))

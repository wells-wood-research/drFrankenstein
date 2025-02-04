import numpy as np
import psi4
from rdkit import Chem


def run_qm_dih_PES_withopt(
    xyz,
    dihedral_atom_index=[1, 2, 3, 4],
    basis="scf/cc-pcvdz",
    FIX_DIH_OPTION="FIXED_DIHEDRAL",
    n_threads=4,
):
    qmol = psi4.core.Molecule.from_string(xyz)
    geom_string = qmol.create_psi4_string_from_molecule()
    ch2 = psi4.geometry(geom_string)
    psi4.set_options(
        {"reference": "rks", "intrafrag_step_limit": 0.1, "freeze_core": "true"}
    )
    if FIX_DIH_OPTION != "FROZEN_DIHEDRAL":
        upper_bound = np.arange(0, 360, 10)
        lower_bound = np.arange(-0.01, 359.99, 10)
    else:
        upper_bound, lower_bound = [0], [0]  # only run one iteration

    PES = []

    psi4.core.set_local_option("OPTKING", "opt_coordinates", "cartesian")
    psi4.core.set_output_file("output_qm.dat", True)
    psi4.set_module_options(
        "optking", {"g_convergence": "gau_loose", "geom_maxiter": 100}
    )
    psi4.core.set_num_threads(nthread=n_threads, quiet=True)
    for lower, upper in zip(lower_bound, upper_bound):
        print(f"Scanning {upper}")

        match FIX_DIH_OPTION:
            case "FROZEN_DIHEDRAL":
                dihedral_string = f"{' '.join(str(i) for i in dihedral_atom_index)}"

                psi4.core.set_local_option(
                    "optking", "FROZEN_DIHEDRAL", dihedral_string
                )
            case "RANGED_DIHEDRAL":
                dihedral_string_ranged = (
                    f"{' '.join(str(i) for i in dihedral_atom_index)} {lower} {upper}"
                )
                psi4.core.set_local_option(
                    "optking", "RANGED_DIHEDRAL", dihedral_string_ranged
                )
            case "FIXED_DIHEDRAL":
                dihedral_string = (
                    f"{' '.join(str(i) for i in dihedral_atom_index)}  {upper}"
                )
                print(dihedral_string)
                psi4.core.set_local_option("optking", "FIXED_DIHEDRAL", dihedral_string)

        try:
            E = psi4.optimize(basis, molecule=ch2) * psi4.constants.hartree2kcalmol
            print(f"Energy on failure : {E}")
            PES.append((upper, E))
            ch2.update_geometry()
        except:
            pass

        PES.append((upper, E, ch2.geometry().np))
        print(PES)
    return PES


from rdkit.Chem import AllChem


def setup_all_qm(mol, dihedrals, n_conf_per_dih):
    AllChem.EmbedMultipleConfs(mol, numConfs=n_conf_per_dih)
    AllChem.MMFFOptimizeMolecule(mol)
    for dihedral_indexes in dihedrals:
        for conf_id in range(n_conf_per_dih):
            conf = mol.GetConformer(conf_id)
            string_xyz = ""
            for atom in mol.GetAtoms():
                atom_idx = atom.GetIdx()
                string_xyz += f" { atom.GetSymbol()} { ' '.join([str(i) for i in conf.GetAtomPosition(atom_idx)])} \n"

            run_qm_dih_PES_withopt(
                string_xyz, dihedral_atom_index=dihedral_indexes, optimise=True
            )


def run_QM_single_point(xyz, basis="scf/6-31G", n_threads=4, optimise=False):
    qmol = psi4.core.Molecule.from_string(xyz)
    geom_string = qmol.create_psi4_string_from_molecule()
    ch2 = psi4.geometry(geom_string)
    psi4.set_options(
        {"reference": "rks", "intrafrag_step_limit": 0.1, "freeze_core": "true"}
    )

    psi4.core.set_local_option("optking", "opt_coordinates", "cartesian")
    psi4.core.set_output_file("output_qm.dat", True)
    psi4.core.set_local_option(
        "optking", {"g_convergence": "gau_loose", "geom_maxiter": 100}
    )
    psi4.core.set_num_threads(nthread=n_threads, quiet=True)
    if optimise == True:
        E = psi4.optimize(basis, molecule=ch2) * psi4.constants.hartree2kcalmol
        print(f"Energy on failure : {E}")
        ch2.update_geometry()
    else:
        E = psi4.energy(basis, molecule=ch2) * psi4.constants.hartree2kcalmol
    return (E, ch2.geometry().np)


def run_some_conf(mol, n_confs, optimise=False):
    AllChem.EmbedMultipleConfs(mol, numConfs=n_confs)
    AllChem.MMFFOptimizeMolecule(mol)
    data = []
    for conf_id in range(n_confs):
        conf = mol.GetConformer(conf_id)
        string_xyz = ""
        for atom in mol.GetAtoms():
            atom_idx = atom.GetIdx()
            string_xyz += f" { atom.GetSymbol()} { ' '.join([str(i) for i in conf.GetAtomPosition(atom_idx)])} \n"
        if conf_id > 1:
            print(string_xyz == test)
        test = string_xyz

        data.append(run_QM_single_point(string_xyz, optimise=optimise))
    return data


def get_xyz_fromrdkitmol(mol):
    xyz = ""
    for i, atom in enumerate(mol.GetAtoms()):
        # positions = molecule.GetConformer().GetAtomPosition(0)
        positions = mol.GetConformer().GetAtomPosition(i)
        xyz += f"{atom.GetSymbol()} {positions.x} {positions.y} {positions.z} \n"
    return xyz


# def find_all_dihedral(mol , frag):


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    from prep import build_molecule_from_smiles

    smile = "CNC(=O)[C@H](CCC[C@@H](NC(C)=O)C(=O)NC)NC(C)=O"
    mol, res_type, backbone_list, capping_list = build_molecule_from_smiles(smile)
    xyz = get_xyz_fromrdkitmol(mol)
    print(backbone_list)

    PES = run_qm_dih_PES_withopt(xyz, basis="HF/6-31G")

    angle = [val[0] for val in PES]
    energy = [val[1] for val in PES]

    plt.plot(angle, energy)

    import pickle as pkl

    pkl.dump(PES, open("example_pes.pkl", "wb"))

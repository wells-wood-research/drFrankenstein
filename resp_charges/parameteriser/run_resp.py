from dataclasses import dataclass
from typing import Literal, Optional

import psi4
import psiresp


@dataclass
class ConformerOptions:
    keep_original_conformer: bool
    minimize: Optional[Literal["uff", "mmff94"]]
    n_conformer_pool: int
    n_max_conformers: int = 200
    rms_tolerance: float = 0.05
    energy_window: float = 30


@dataclass
class RESP_Params:
    VDW_SCALE_FACTORS = [1.4, 1.6, 1.8, 2.0]
    VDW_POINT_DENSITY = 1.0
    RESP_A = 0.0005
    RESP_B = 0.1
    BASIS_ESP = "6-31G**"  # "HF-CEP-31G"


def run_resp(geometry, parameters: RESP_Params):
    import resp

    mol = psi4.geometry(geometry)
    mol.update_geometry()

    options = {
        "VDW_SCALE_FACTORS": parameters.VDW_SCALE_FACTORS,
        "VDW_POINT_DENSITY": 1.0,
        "RESP_A": 0.0005,
        "RESP_B": 0.1,
        "BASIS_ESP": parameters.BASIS_ESP,
    }

    # Call for first stage fit
    charges1 = resp.resp([mol], options)
    print("Electrostatic Potential Charges")
    print(charges1[0])
    print("Restrained Electrostatic Potential Charges")
    print(charges1[1])

    # Change the value of the RESP parameter A
    options["RESP_A"] = 0.001

    # Add constraint for atoms fixed in second stage fit
    constraint_charge = []
    for i in range(4, 8):
        constraint_charge.append([charges1[1][i], [i + 1]])
    options["constraint_charge"] = constraint_charge
    options["constraint_group"] = [[2, 3, 4]]
    options["grid"] = ["1_%s_grid.dat" % mol.name()]
    options["esp"] = ["1_%s_grid_esp.dat" % mol.name()]

    # Call for second stage fit
    charges2 = resp.resp([mol], options)

    # Get RESP charges
    print("\nStage Two:\n")
    print("RESP Charges")
    print(charges2[1])


def impose_amber14ff_charges(list_atoms):
    return {}



if __name__ == "__main__":
    from prep import build_molecule_from_smiles

    smile = "CNC(=O)[C@H](CCC[C@@H](NC(C)=O)C(=O)NC)NC(C)=O"
    mol, _, backbone_list, capping_list = build_molecule_from_smiles(smile)
    print(capping_list)
    constraint_capping = run_psiresp(mol, backbone_list, capping_list)
    without_constraint = run_psiresp(mol, backbone_list, [])
    print(constraint_capping)
    atom_names = [
        f"{x}{i}"
        for i, x in enumerate([atom.GetSymbol() for atom in mol.GetAtoms()], 1)
    ]
    import pandas as pd

    df = pd.DataFrame(
        {
            "Atom names": atom_names,
            "RESP": constraint_capping,
            "without constraint": without_constraint,
        }
    )

    df.set_index("Atom names", inplace=True)

    df.plot(
        ylabel="Charge",
        title="with and without restraint",
        kind="bar",
    )
    import matplotlib.pyplot as plt

    plt.show()

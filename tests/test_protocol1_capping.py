import unittest

import numpy as np
import pandas as pd

from Laboratory.Experiments.Protocol_1_Capping import Capping_Assistant
from Laboratory.Experiments.Protocol_1_Capping import Capping_Builders
from Laboratory.Experiments.Protocol_1_Capping import Capping_Geometry as geom


class TestProtocol1Capping(unittest.TestCase):
    def test_n_terminus_deletes_enough_hydrogens_for_trivalent_n(self):
        config = {
            "moleculeInfo": {
                "backboneAliases": {
                    "H": ["H"],
                }
            }
        }
        bonded_atoms = ["CA", "H", "H2", "H3"]

        atoms_to_delete = Capping_Assistant.decide_atom_to_delete_N_termini(bonded_atoms, config)

        self.assertEqual(set(atoms_to_delete), {"H2", "H3"})

    def test_ace_carbonyl_improper_is_trans_180(self):
        mol_df = pd.DataFrame(
            {
                "ATOM_ID": [1, 2, 3, 4, 5],
                "ATOM_NAME": ["N", "H", "CA", "C", "O"],
                "RES_NAME": ["ALA", "ALA", "ALA", "ALA", "ALA"],
                "RES_ID": [1, 1, 1, 1, 1],
                "X": [0.0, -0.5, 1.5, 2.5, 3.0],
                "Y": [0.0, -0.5, 0.0, 1.0, 2.0],
                "Z": [0.0, 0.5, 0.0, 0.0, 0.0],
                "OCCUPANCY": [1.0] * 5,
                "TEMP_FACTOR": [0.0] * 5,
                "ELEMENT": ["N", "H", "C", "C", "O"],
            }
        )
        ace_template_df = pd.DataFrame(
            {
                "ATOM_NAME": ["C_C", "O_C", "C2_C", "H1_C", "H2_C", "H3_C"],
                "X": [0.0, 1.2, -1.5, -2.0, -1.2, -1.8],
                "Y": [0.0, 0.0, 0.0, 0.9, -0.8, -0.1],
                "Z": [0.0, 0.0, 0.0, 0.0, 0.5, -0.7],
            }
        )

        ace_cap = Capping_Builders.build_ace_cap(mol_df, "N", ace_template_df)
        validation = Capping_Builders.validate_geometry(mol_df, ace_cap, "N", "N")

        self.assertAlmostEqual(validation["N-C_C-O_C_angle"], 120.0, places=5)
        self.assertAlmostEqual(validation["N-C_C-C2_C_angle"], 120.0, places=5)
        self.assertAlmostEqual(validation["O_C-C_C-C2_C_angle"], 120.0, places=5)
        self.assertIn("O_C-C_C-N-C2_C_dihedral", validation)
        self.assertAlmostEqual(abs(validation["O_C-C_C-N-C2_C_dihedral"]), 180.0, places=5)
        self.assertIn("O_C-C_C-N-CA_dihedral", validation)
        self.assertAlmostEqual(abs(validation["O_C-C_C-N-CA_dihedral"]), 180.0, places=5)

    def test_nme_nitrogen_angles_are_trigonal_planar(self):
        mol_df = pd.DataFrame(
            {
                "ATOM_ID": [1, 2, 3, 4, 5],
                "ATOM_NAME": ["N", "H", "CA", "C", "O"],
                "RES_NAME": ["ALA", "ALA", "ALA", "ALA", "ALA"],
                "RES_ID": [1, 1, 1, 1, 1],
                "X": [0.0, -0.5, 1.5, 2.5, 3.0],
                "Y": [0.0, -0.5, 0.0, 1.0, 2.0],
                "Z": [0.0, 0.5, 0.0, 0.0, 0.0],
                "OCCUPANCY": [1.0] * 5,
                "TEMP_FACTOR": [0.0] * 5,
                "ELEMENT": ["N", "H", "C", "C", "O"],
            }
        )
        nme_template_df = pd.DataFrame(
            {
                "ATOM_NAME": ["N_N", "H_N", "C_N", "H1_N", "H2_N", "H3_N"],
                "X": [0.0, 0.0, 1.4, 1.9, 1.6, 1.8],
                "Y": [0.0, 1.0, 0.0, 0.8, -0.9, 0.1],
                "Z": [0.0, 0.0, 0.0, -0.6, 0.4, 0.9],
            }
        )

        nme_cap = Capping_Builders.build_nme_cap(mol_df, "C", nme_template_df)
        validation = Capping_Builders.validate_geometry(mol_df, nme_cap, "C", "C")

        self.assertAlmostEqual(validation["C-N_N-H_N_angle"], 120.0, places=5)
        self.assertAlmostEqual(validation["C-N_N-C_N_angle"], 120.0, places=5)
        self.assertAlmostEqual(validation["H_N-N_N-C_N_angle"], 120.0, places=5)

    def test_cap_clash_resolver_rotates_rigidly_around_attachment_axis(self):
        mol_df = pd.DataFrame(
            {
                "ATOM_NAME": ["T", "X"],
                "X": [0.0, 0.0],
                "Y": [0.0, 2.0],
                "Z": [0.0, 0.0],
            }
        )
        cap_df = pd.DataFrame(
            {
                "ATOM_NAME": ["P", "A", "B"],
                "X": [1.0, 1.0, 1.0],
                "Y": [0.0, 2.0, 0.0],
                "Z": [0.0, 0.0, 1.0],
            }
        )

        before_angle = geom.calculate_angle(
            geom.get_coords(mol_df, "T"),
            geom.get_coords(cap_df, "P"),
            geom.get_coords(cap_df, "A"),
        )
        before_tp = geom.calculate_distance(geom.get_coords(mol_df, "T"), geom.get_coords(cap_df, "P"))

        resolved = Capping_Builders.resolve_cap_nonbonded_clashes(
            mol_df=mol_df,
            cap_df=cap_df,
            terminal_atom="T",
            pivot_atom="P",
            min_distance=1.5,
            rotation_steps=72,
        )

        blocker = geom.get_coords(mol_df, "X")
        resolved_coords = resolved[["X", "Y", "Z"]].astype(float).values
        min_distance = np.min(np.linalg.norm(resolved_coords - blocker[None, :], axis=1))

        after_angle = geom.calculate_angle(
            geom.get_coords(mol_df, "T"),
            geom.get_coords(resolved, "P"),
            geom.get_coords(resolved, "A"),
        )
        after_tp = geom.calculate_distance(geom.get_coords(mol_df, "T"), geom.get_coords(resolved, "P"))

        self.assertGreaterEqual(min_distance, 1.5 - 1e-6)
        self.assertAlmostEqual(before_angle, after_angle, places=5)
        self.assertAlmostEqual(before_tp, after_tp, places=5)

    def test_ace_builder_avoids_nonbonded_cap_clashes(self):
        mol_df = pd.DataFrame(
            {
                "ATOM_ID": [1, 2, 3, 4, 5],
                "ATOM_NAME": ["N", "H", "CA", "C", "O"],
                "RES_NAME": ["ALA", "ALA", "ALA", "ALA", "ALA"],
                "RES_ID": [1, 1, 1, 1, 1],
                "X": [0.0, -0.5, 1.5, 2.5, 3.0],
                "Y": [0.0, -0.5, 0.0, 1.0, 2.0],
                "Z": [0.0, 0.5, 0.0, 0.0, 0.0],
                "OCCUPANCY": [1.0] * 5,
                "TEMP_FACTOR": [0.0] * 5,
                "ELEMENT": ["N", "H", "C", "C", "O"],
            }
        )
        ace_template_df = pd.DataFrame(
            {
                "ATOM_NAME": ["C_C", "O_C", "C2_C", "H1_C", "H2_C", "H3_C"],
                "X": [0.0, 1.2, -1.5, -2.0, -1.2, -1.8],
                "Y": [0.0, 0.0, 0.0, 0.9, -0.8, -0.1],
                "Z": [0.0, 0.0, 0.0, 0.0, 0.5, -0.7],
            }
        )

        ace_without_blocker = Capping_Builders.build_ace_cap(mol_df, "N", ace_template_df)
        blocker = geom.get_coords(ace_without_blocker, "H1_C")

        blocker_row = pd.DataFrame(
            {
                "ATOM_ID": [99],
                "ATOM_NAME": ["X"],
                "RES_NAME": ["ALA"],
                "RES_ID": [1],
                "X": [float(blocker[0])],
                "Y": [float(blocker[1])],
                "Z": [float(blocker[2])],
                "OCCUPANCY": [1.0],
                "TEMP_FACTOR": [0.0],
                "ELEMENT": ["C"],
            }
        )
        mol_blocked = pd.concat([mol_df, blocker_row], ignore_index=True)

        ace_with_blocker = Capping_Builders.build_ace_cap(mol_blocked, "N", ace_template_df)
        validation = Capping_Builders.validate_geometry(mol_blocked, ace_with_blocker, "N", "N")

        self.assertGreaterEqual(validation["min_cap_nonbonded_distance"], 1.5 - 1e-6)
        self.assertAlmostEqual(validation["N-C_C-O_C_angle"], 120.0, places=5)
        self.assertAlmostEqual(validation["N-C_C-C2_C_angle"], 120.0, places=5)
        self.assertAlmostEqual(validation["O_C-C_C-C2_C_angle"], 120.0, places=5)


if __name__ == "__main__":
    unittest.main()

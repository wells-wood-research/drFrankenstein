from pathlib import Path

import parmed

# import sire as sr
from psiomm import molecule
from psiomm import psi_omm as po


def get_decomposed_energy_sire(topology: Path, trajectory: Path) -> dict:
    mols = sr.load(topology, "SYSTEM.crd")
    energy_decomp = mols[0].energy().components()

    return energy_decomp


def get_decomposed_energy_parmed(topology: Path, coordinates):
    parmed_top = parmed.load(topology)
    context = create_context(parmed_top, coordinates)
    energy_decomp = []
    for coord in coordinates:
        context.setPositions(coord)
        energy_decomp.append(parmed.openmm.energy_decomposition(parmed_top, context))

    return energy_decomp


def create_context(topology: parmed.Topology, coordinates):
    system = topology.createSystem()
    return system.context


def energy_without_dih(energy_decomp: list):
    energies = []
    for e_d in energy_decomp:
        energies.append(0)
        for key in e_d.keys():
            if key != "DIHEDRAL_TERM":
                energies[-1] += e_d[key]
    return energies


def tot_energy(energy_decomp: list):
    energies = []
    for e_d in energy_decomp:
        energies.append(0)
        for key in e_d.keys():
            energies[-1] += e_d[key]
    return energies


def get_topology(xyzstring):
    ### Promising need to test if topology compatible with parmed
    solute = molecule.Molecule.from_xyz_string(xyzstring)
    solute.generate_bonds()
    solute.generate_atom_types()
    solute.generate_charges()

    # check generated atoms and charges
    for i, at in enumerate(solute.atom_types):
        print(i, at, solute.charges[i])

    topology = po.make_topology(solute)
    forcefield = ForceField("gaff2.xml", "tip3p.xml")

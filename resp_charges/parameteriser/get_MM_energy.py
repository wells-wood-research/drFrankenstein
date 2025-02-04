import openmm as mm
import parmed
from openmm import app


def return_MM(topology_file, frames):
    # Instantiate the parm and create the system

    parm = parmed.load_file(topology_file)

    system = parm.createSystem(
        nonbondedMethod=app.PME, nonbondedCutoff=8 * pmd.unit.angstrom
    )

    # Make the context and set the positions
    context = mm.Context(system, mm.VerletIntegrator(0.001))

    energy_total = []
    for positions in frames:
        context.setPositions(positions)
        # Find the energy decomposition
        energy_total.append(parmed.openmm.energy_decomposition(parm, context)["total"])
    return energy_total

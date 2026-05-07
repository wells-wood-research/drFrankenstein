import numpy as np


## CLEAN CODE ##
class FilePath:
    pass


def get_ordered_conformer_xyzs(config: dict, temperature: float = 300) -> list[FilePath]:
    """
    Returns all conformers in Boltzmann-weighted random order (no replacement).
    """
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    conformerEnergies = config["runtimeInfo"]["madeByConformers"].get("conformerEnergies")
    if len(conformerXyzs) <= 1:
        return conformerXyzs
    if not conformerEnergies or len(conformerEnergies) != len(conformerXyzs):
        return conformerXyzs

    seed = config["miscInfo"]["seed"]
    rng = np.random.default_rng(seed)
    kBoltzmann = 0.0019872041  # Boltzmann constant in kcal/mol·K
    kT = kBoltzmann * temperature

    sortedKeys = sorted(conformerEnergies.keys())
    sortedConformerXyzs = [conformerXyzs[list(conformerEnergies.keys()).index(key)] for key in sortedKeys]
    energies = [conformerEnergies[key] for key in sortedKeys]
    minEnergy = min(energies)
    shiftedEnergies = [energy - minEnergy for energy in energies]

    boltzmannFactors = [np.exp(-energy / kT) for energy in shiftedEnergies]
    totalWeight = sum(boltzmannFactors)
    probabilities = [factor / totalWeight for factor in boltzmannFactors]

    selectedIndices = rng.choice(
        len(sortedConformerXyzs),
        size=len(sortedConformerXyzs),
        replace=False,
        p=probabilities
    )

    return [sortedConformerXyzs[i] for i in selectedIndices]


def select_conformer_xyzs(config: dict, nConformers: int, temperature: float = 300) -> list[FilePath]:
    """
    Returns conformers selected from a seeded Boltzmann-weighted ordering.
    """
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    if nConformers == -1 or nConformers >= len(conformerXyzs):
        return conformerXyzs

    orderedConformerXyzs = get_ordered_conformer_xyzs(config, temperature=temperature)
    return orderedConformerXyzs[:nConformers]

import numpy as np

## CLEAN CODE ##
class FilePath:
    pass


_ENERGY_METHOD = "ENERGY"
_DIVERSE_METHOD = "DIVERSE"
_LEGACY_BOLTZMANN_METHOD = "BOLTZMANN"
_LEGACY_PCA_KMEANS_METHOD = "PCA_KMEANS"
_LEGACY_PCA_KNN_METHOD = "PCA_KNN"


def _normalize_selection_method(selection_method: str) -> str:
    if not selection_method:
        return _ENERGY_METHOD
    normalized = selection_method.upper()
    if normalized in (_ENERGY_METHOD, _DIVERSE_METHOD):
        return normalized
    if normalized == _LEGACY_BOLTZMANN_METHOD:
        return _ENERGY_METHOD
    if normalized in (_LEGACY_PCA_KMEANS_METHOD, _LEGACY_PCA_KNN_METHOD):
        return _DIVERSE_METHOD
    return selection_method


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
    Returns conformers selected using the miscInfo.conformerSelectionMethods strategy.
    """
    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    if nConformers == -1 or nConformers >= len(conformerXyzs):
        return conformerXyzs

    selection_method = config.get("miscInfo", {}).get("conformerSelectionMethods")
    if selection_method is None:
        selection_method = config.get("chargeFittingInfo", {}).get("conformerSelectionMethod", _ENERGY_METHOD)
    selection_method = _normalize_selection_method(selection_method)

    if selection_method == _ENERGY_METHOD:
        orderedConformerXyzs = get_ordered_conformer_xyzs(config, temperature=temperature)
        return orderedConformerXyzs[:nConformers]
    if selection_method == _DIVERSE_METHOD:
        return _select_conformers_by_pca_kmeans(config, nConformers)

    raise ValueError(
        f"Unknown conformer selection method '{selection_method}'. "
        f"Expected '{_ENERGY_METHOD}' or '{_DIVERSE_METHOD}'."
    )


def _select_conformers_by_pca_kmeans(config: dict, nConformers: int) -> list[FilePath]:
    from os import path as p
    import pandas as pd
    from sklearn.cluster import KMeans

    from Experiments.Protocol_8_Reporter import Conformer_PCA_Monster

    conformerXyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    if nConformers <= 0:
        return []
    if nConformers >= len(conformerXyzs):
        return conformerXyzs

    conformerEnergies = config["runtimeInfo"]["madeByConformers"].get("conformerEnergies")
    if not conformerEnergies or len(conformerEnergies) != len(conformerXyzs):
        raise ValueError("PCA selection requires conformer energies for every conformer.")

    twisting_info = config.get("runtimeInfo", {}).get("madeByTwisting", {})
    if "allDihedrals" not in twisting_info:
        raise ValueError("PCA selection requires runtimeInfo.madeByTwisting.allDihedrals.")

    torsions = Conformer_PCA_Monster._flatten_all_dihedrals(config)
    if not torsions:
        raise ValueError("PCA selection requires torsion definitions from allDihedrals.")

    angle_rows = [Conformer_PCA_Monster._measure_torsion_angles(conformer_xyz, torsions) for conformer_xyz in conformerXyzs]
    angle_df = pd.DataFrame(angle_rows)
    pca_df, _, _ = Conformer_PCA_Monster._run_pca(angle_df)
    pca_df["energy"] = [
        Conformer_PCA_Monster._lookup_conformer_energy(conformer_name, conformerEnergies)
        for conformer_name in pca_df["conformer"]
    ]

    rng_seed = config.get("miscInfo", {}).get("seed", 1818)
    kmeans = KMeans(n_clusters=nConformers, random_state=rng_seed, n_init="auto")
    pca_coords = pca_df[["pc1", "pc2"]].to_numpy()
    pca_df["cluster"] = kmeans.fit_predict(pca_coords)

    selected_names = []
    for cluster_id in range(nConformers):
        cluster_rows = pca_df[pca_df["cluster"] == cluster_id]
        if cluster_rows.empty:
            continue
        min_row = cluster_rows.loc[cluster_rows["energy"].astype(float).idxmin()]
        selected_names.append(min_row["conformer"])

    if len(selected_names) < nConformers:
        remaining = pca_df[~pca_df["conformer"].isin(selected_names)].sort_values("energy")
        for conformer_name in remaining["conformer"].tolist():
            if len(selected_names) >= nConformers:
                break
            selected_names.append(conformer_name)

    name_to_xyz = {p.splitext(p.basename(conformer_xyz))[0]: conformer_xyz for conformer_xyz in conformerXyzs}
    return [name_to_xyz[conformer_name] for conformer_name in selected_names if conformer_name in name_to_xyz]

import os
from os import path as p

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from . import Reporting_Assistant
from ..Protocol_5_Twisting import Twisted_Assistant


def _flatten_all_dihedrals(config: dict) -> list[dict]:
    """Flatten the twisted dihedral structure into a list of torsion records."""
    all_dihedrals = config["runtimeInfo"]["madeByTwisting"]["allDihedrals"]
    torsions = []

    for torsion_group, torsion_map in all_dihedrals.items():
        for torsion_tag, torsion_data in torsion_map.items():
            torsions.append(
                {
                    "torsionGroup": torsion_group,
                    "torsionTag": torsion_tag,
                    "torsionLabel": f"{torsion_group}:{torsion_tag}",
                    "atomIndexes": torsion_data["ATOM_INDEXES"],
                }
            )

    return torsions


def _measure_torsion_angles(conformer_xyz: str, torsions: list[dict]) -> dict:
    """Measure every tracked torsion angle for a conformer."""
    atom_coords = Twisted_Assistant.xyz_to_np_array(conformer_xyz)
    conformer_name = p.splitext(p.basename(conformer_xyz))[0]
    measured_angles = {"conformer": conformer_name}

    for torsion in torsions:
        angle = Twisted_Assistant.calculate_torsion_angle(atom_coords, torsion["atomIndexes"])
        measured_angles[torsion["torsionLabel"]] = angle

    return measured_angles


def _lookup_conformer_energy(conformer_name: str, conformer_energies: dict) -> float:
    """Resolve a conformer's energy from the energy lookup table."""
    if conformer_name in conformer_energies:
        return float(conformer_energies[conformer_name])

    suffix = p.basename(conformer_name)
    if suffix in conformer_energies:
        return float(conformer_energies[suffix])

    raise KeyError(f"Energy not found for conformer '{conformer_name}'.")


def _build_conformer_html_map(conformer_names: list[str], images_dir: str, reporter_dir: str) -> dict:
    """Map conformer names to their generated HTML files."""
    conformer_images_dir = p.join(images_dir, "conformer_images")
    conformer_html_map = {}
    for conformer_name in conformer_names:
        html_path = p.join(conformer_images_dir, f"{conformer_name}.html")
        if p.exists(html_path):
            conformer_html_map[conformer_name] = p.relpath(html_path, reporter_dir)
    return conformer_html_map


def _conformer_names_from_xyzs(conformer_xyzs: list[str]) -> list[str]:
    """Extract conformer names from XYZ file paths."""
    return [p.splitext(p.basename(conformer_xyz))[0] for conformer_xyz in conformer_xyzs]


def _infer_torsion_scan_xyzs(config: dict, conformer_xyzs: list[str]) -> list[str]:
    """Infer which conformers were used for torsion scans."""
    n_conformers_requested = config.get("torsionScanInfo", {}).get("nConformers", -1)
    if n_conformers_requested == -1 or n_conformers_requested >= len(conformer_xyzs):
        return conformer_xyzs

    from OperatingTools import select_conformers

    ordered_conformers = select_conformers.get_ordered_conformer_xyzs(config)
    return ordered_conformers[:n_conformers_requested]


def _run_pca(angle_df: pd.DataFrame) -> tuple[pd.DataFrame, list[float], bool]:
    """Run PCA on torsion angle data and return the transformed coordinates."""
    feature_columns = [column for column in angle_df.columns if column != "conformer"]
    feature_matrix = angle_df[feature_columns].to_numpy()

    if feature_matrix.shape[0] < 2 or feature_matrix.shape[1] == 0:
        pca_df = pd.DataFrame(
            {
                "conformer": angle_df["conformer"],
                "pc1": np.arange(len(angle_df), dtype=float),
                "pc2": np.zeros(len(angle_df), dtype=float),
            }
        )
        return pca_df, [1.0, 0.0], False

    scaled_features = StandardScaler().fit_transform(feature_matrix)
    n_components = min(2, scaled_features.shape[0], scaled_features.shape[1])
    pca = PCA(n_components=n_components)
    transformed = pca.fit_transform(scaled_features)

    if n_components == 1:
        transformed = np.column_stack([transformed[:, 0], np.zeros(transformed.shape[0])])
        explained_variance = [float(pca.explained_variance_ratio_[0]), 0.0]
    else:
        explained_variance = [float(value) for value in pca.explained_variance_ratio_[:2]]

    pca_df = pd.DataFrame(
        {
            "conformer": angle_df["conformer"],
            "pc1": transformed[:, 0],
            "pc2": transformed[:, 1],
        }
    )
    return pca_df, explained_variance, True


def _build_pca_plot_html(
    pca_df: pd.DataFrame,
    explained_variance: list[float],
    selected_for_charges: set[str] | None = None,
    selected_for_scans: set[str] | None = None,
) -> str:
    """Build the PCA scatter plot HTML."""
    energy_values = pca_df["energy"].astype(float).tolist()
    colorscale, cmin, cmax = _build_energy_colorscale(energy_values)
    selected_for_charges = selected_for_charges or set()
    selected_for_scans = selected_for_scans or set()

    fig = go.Figure(
        data=[
            go.Scatter(
                x=pca_df["pc1"],
                y=pca_df["pc2"],
                mode="markers",
                showlegend=False,
                customdata=np.column_stack((pca_df["conformer"], pca_df["energy"])),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "PC1: %{x:.4f}<br>"
                    "PC2: %{y:.4f}<br>"
                    "Energy: %{customdata[1]:.4f} kcal/mol<extra></extra>"
                ),
                marker=dict(
                    size=14,
                    color=pca_df["energy"],
                    colorscale=colorscale,
                    cmin=cmin,
                    cmax=cmax,
                    colorbar=dict(
                        title=dict(text="Relative Energy\n(kcal/mol)", font=dict(color="yellow")),
                        tickfont=dict(color="yellow"),
                        outlinecolor="yellow",
                        x=1.05,
                        y=0.4,
                        len=0.8,
                        thickness=16,
                    ),
                    line=dict(color="black", width=1),
                ),
            )
        ]
    )

    scan_selected_df = pca_df[pca_df["conformer"].isin(selected_for_scans)]
    if not scan_selected_df.empty:
        fig.add_trace(
            go.Scatter(
                x=scan_selected_df["pc1"],
                y=scan_selected_df["pc2"],
                mode="markers",
                name="Torsion scan",
                customdata=np.column_stack((scan_selected_df["conformer"], scan_selected_df["energy"])),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "PC1: %{x:.4f}<br>"
                    "PC2: %{y:.4f}<br>"
                    "Energy: %{customdata[1]:.4f} kcal/mol<extra></extra>"
                ),
                marker=dict(
                    size=18,
                    color="#00FFFF",
                    symbol="square-open",
                    line=dict(color="#00FFFF", width=2),
                ),
            )
        )

    charge_selected_df = pca_df[pca_df["conformer"].isin(selected_for_charges)]
    if not charge_selected_df.empty:
        fig.add_trace(
            go.Scatter(
                x=charge_selected_df["pc1"],
                y=charge_selected_df["pc2"],
                mode="markers",
                name="Charge fitting",
                customdata=np.column_stack((charge_selected_df["conformer"], charge_selected_df["energy"])),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "PC1: %{x:.4f}<br>"
                    "PC2: %{y:.4f}<br>"
                    "Energy: %{customdata[1]:.4f} kcal/mol<extra></extra>"
                ),
                marker=dict(
                    size=18,
                    color="#FFA500",
                    symbol="diamond-open",
                    line=dict(color="#FFA500", width=2),
                ),
            )
        )

    fig.update_layout(
        title=dict(text="Conformer PCA of torsional space", font=dict(color="yellow")),
        template="plotly_dark",
        paper_bgcolor="#111111",
        plot_bgcolor="#111111",
        font=dict(color="yellow", family="Consolas, Courier New, monospace"),
        margin=dict(l=70, r=150, t=60, b=60),
        legend=dict(
            x=1.05,
            y=0.98,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(0, 0, 0, 0.6)",
            bordercolor="yellow",
            borderwidth=1,
            font=dict(color="yellow"),
        ),
    )
    fig.update_xaxes(
        title=dict(text=f"PC1 ({explained_variance[0] * 100:.1f}% variance)", font=dict(color="yellow")),
        tickfont=dict(color="yellow"),
        color="yellow",
        gridcolor="#444444",
        zerolinecolor="#888888",
        linecolor="yellow",
    )
    fig.update_yaxes(
        title=dict(text=f"PC2 ({explained_variance[1] * 100:.1f}% variance)", font=dict(color="yellow")),
        tickfont=dict(color="yellow"),
        color="yellow",
        gridcolor="#444444",
        zerolinecolor="#888888",
        linecolor="yellow",
    )

    return fig.to_html(full_html=False, include_plotlyjs="cdn", config={"displayModeBar": False})


def _build_energy_colorscale(energy_values: list[float]) -> tuple[list[list[float | str]], float, float]:
    """Build a Plotly colorscale for energy values."""
    vibrant_colors = Reporting_Assistant.make_vibrant_colors()
    green = vibrant_colors[1]
    magenta = vibrant_colors[4]
    colorscale = [
        [0.0, green],
        [0.5, "#FFFFFF"],
        [1.0, magenta],
    ]

    cmin = min(energy_values)
    cmax = max(energy_values)
    if cmin == cmax:
        cmax = cmin + 1.0

    return colorscale, cmin, cmax


def _backbone_aliases_ready(backbone_aliases: dict | None) -> bool:
    """Check whether the backbone aliases needed for Ramachandran plots are present."""
    if not isinstance(backbone_aliases, dict):
        return False
    required = ("C", "N", "O", "CA")
    return all(backbone_aliases.get(key) for key in required)


def _classify_phi_psi_dihedrals(phi_psi_dihedrals: dict, backbone_aliases: dict) -> tuple[dict[int, list[int]], dict[int, list[int]]]:
    """Classify phi and psi dihedrals by CA index."""
    n_aliases = set(backbone_aliases.get("N", []))
    ca_aliases = set(backbone_aliases.get("CA", []))
    c_aliases = set(backbone_aliases.get("C", []))

    phi_by_ca = {}
    psi_by_ca = {}

    for dihedral_data in phi_psi_dihedrals.values():
        atom_names = dihedral_data.get("ATOM_NAMES", [])
        atom_indexes = dihedral_data.get("ATOM_INDEXES", [])
        if len(atom_names) < 4 or len(atom_indexes) < 4:
            continue
        center_a = atom_names[1]
        center_b = atom_names[2]

        if (center_a in n_aliases and center_b in ca_aliases) or (center_b in n_aliases and center_a in ca_aliases):
            ca_index = atom_indexes[2] if center_b in ca_aliases else atom_indexes[1]
            phi_by_ca[ca_index] = atom_indexes
            continue
        if (center_a in ca_aliases and center_b in c_aliases) or (center_b in ca_aliases and center_a in c_aliases):
            ca_index = atom_indexes[1] if center_a in ca_aliases else atom_indexes[2]
            psi_by_ca[ca_index] = atom_indexes

    return phi_by_ca, psi_by_ca


def _measure_ramachandran_angles(conformer_xyzs: list[str], phi_psi_pairs: list[tuple[int, list[int], list[int]]], conformer_energies: dict) -> list[dict]:
    """Measure phi/psi angles for each conformer and residue pair."""
    ramachan_rows = []
    for conformer_xyz in conformer_xyzs:
        atom_coords = Twisted_Assistant.xyz_to_np_array(conformer_xyz)
        conformer_name = p.splitext(p.basename(conformer_xyz))[0]
        energy = _lookup_conformer_energy(conformer_name, conformer_energies)
        for ca_index, phi_indexes, psi_indexes in phi_psi_pairs:
            phi_angle = Twisted_Assistant.calculate_torsion_angle(atom_coords, phi_indexes)
            psi_angle = Twisted_Assistant.calculate_torsion_angle(atom_coords, psi_indexes)
            ramachan_rows.append(
                {
                    "conformer": conformer_name,
                    "residueKey": f"CA {ca_index}",
                    "phi": phi_angle,
                    "psi": psi_angle,
                    "energy": energy,
                }
            )
    return ramachan_rows


def _build_ramachandran_plot_html(
    rama_df: pd.DataFrame,
    selected_for_charges: set[str] | None = None,
    selected_for_scans: set[str] | None = None,
) -> str:
    """Build the Ramachandran plot HTML."""
    energy_values = rama_df["energy"].astype(float).tolist()
    colorscale, cmin, cmax = _build_energy_colorscale(energy_values)
    selected_for_charges = selected_for_charges or set()
    selected_for_scans = selected_for_scans or set()

    fig = go.Figure(
        data=[
            go.Scatter(
                x=rama_df["phi"],
                y=rama_df["psi"],
                mode="markers",
                showlegend=False,
                customdata=np.column_stack((rama_df["conformer"], rama_df["residueKey"], rama_df["energy"])),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "%{customdata[1]}<br>"
                    "Phi: %{x:.1f}°<br>"
                    "Psi: %{y:.1f}°<br>"
                    "Energy: %{customdata[2]:.4f} kcal/mol<extra></extra>"
                ),
                marker=dict(
                    size=14,
                    color=rama_df["energy"],
                    colorscale=colorscale,
                    cmin=cmin,
                    cmax=cmax,
                    colorbar=dict(
                        title=dict(text="Energy (kcal/mol)", font=dict(color="yellow")),
                        tickfont=dict(color="yellow"),
                        outlinecolor="yellow",
                        x=1.05,
                        y=0.4,
                        len=0.8,
                        thickness=16,
                    ),
                    line=dict(color="black", width=1),
                ),
            )
        ]
    )

    scan_selected_df = rama_df[rama_df["conformer"].isin(selected_for_scans)]
    if not scan_selected_df.empty:
        fig.add_trace(
            go.Scatter(
                x=scan_selected_df["phi"],
                y=scan_selected_df["psi"],
                mode="markers",
                name="Torsion scan",
                customdata=np.column_stack(
                    (scan_selected_df["conformer"], scan_selected_df["residueKey"], scan_selected_df["energy"])
                ),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "%{customdata[1]}<br>"
                    "Phi: %{x:.1f}°<br>"
                    "Psi: %{y:.1f}°<br>"
                    "Energy: %{customdata[2]:.4f} kcal/mol<extra></extra>"
                ),
                marker=dict(
                    size=16,
                    color="#00FFFF",
                    symbol="square-open",
                    line=dict(color="#00FFFF", width=2),
                ),
            )
        )

    charge_selected_df = rama_df[rama_df["conformer"].isin(selected_for_charges)]
    if not charge_selected_df.empty:
        fig.add_trace(
            go.Scatter(
                x=charge_selected_df["phi"],
                y=charge_selected_df["psi"],
                mode="markers",
                name="Charge fitting",
                customdata=np.column_stack(
                    (charge_selected_df["conformer"], charge_selected_df["residueKey"], charge_selected_df["energy"])
                ),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "%{customdata[1]}<br>"
                    "Phi: %{x:.1f}°<br>"
                    "Psi: %{y:.1f}°<br>"
                    "Energy: %{customdata[2]:.4f} kcal/mol<extra></extra>"
                ),
                marker=dict(
                    size=16,
                    color="#FFA500",
                    symbol="diamond-open",
                    line=dict(color="#FFA500", width=2),
                ),
            )
        )

    fig.update_layout(
        title=dict(text="Ramachandran plot (Phi vs Psi)", font=dict(color="yellow")),
        template="plotly_dark",
        paper_bgcolor="#111111",
        plot_bgcolor="#111111",
        font=dict(color="yellow", family="Consolas, Courier New, monospace"),
        margin=dict(l=70, r=150, t=60, b=60),
        legend=dict(
            x=1.05,
            y=0.98,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(0, 0, 0, 0.6)",
            bordercolor="yellow",
            borderwidth=1,
            font=dict(color="yellow"),
        ),
    )
    fig.update_xaxes(
        title=dict(text="Phi (°)", font=dict(color="yellow")),
        tickfont=dict(color="yellow"),
        color="yellow",
        gridcolor="#444444",
        zerolinecolor="#888888",
        linecolor="yellow",
        range=[-180, 180],
        dtick=60,
    )
    fig.update_yaxes(
        title=dict(text="Psi (°)", font=dict(color="yellow")),
        tickfont=dict(color="yellow"),
        color="yellow",
        gridcolor="#444444",
        zerolinecolor="#888888",
        linecolor="yellow",
        range=[-180, 180],
        dtick=60,
    )

    return fig.to_html(full_html=False, include_plotlyjs="cdn", config={"displayModeBar": False})


def process_conformer_pca_results(config: dict) -> dict:
    """Generate PCA and Ramachandran report data for conformers."""
    images_dir = config["runtimeInfo"]["madeByReporting"]["imagesDir"]
    reporter_dir = config["runtimeInfo"]["madeByReporting"]["reporterDir"]
    conformer_dir = p.join(images_dir, "conformer_pca")
    os.makedirs(conformer_dir, exist_ok=True)

    torsions = _flatten_all_dihedrals(config)
    conformer_xyzs = config["runtimeInfo"]["madeByConformers"]["conformerXyzs"]
    conformer_energies = config["runtimeInfo"]["madeByConformers"]["conformerEnergies"]
    if not torsions:
        raise ValueError("No torsions were available in runtimeInfo.madeByTwisting.allDihedrals.")
    if not conformer_xyzs:
        raise ValueError("No conformer XYZ files were available for PCA analysis.")

    angle_rows = [_measure_torsion_angles(conformer_xyz, torsions) for conformer_xyz in conformer_xyzs]
    angle_df = pd.DataFrame(angle_rows)
    pca_df, explained_variance, pca_available = _run_pca(angle_df)
    pca_df["energy"] = [
        _lookup_conformer_energy(conformer_name, conformer_energies)
        for conformer_name in pca_df["conformer"]
    ]
    charges_xyzs = config.get("runtimeInfo", {}).get("madeByCharges", {}).get("conformerXyzsForCharges", [])
    torsion_scan_xyzs = config.get("runtimeInfo", {}).get("madeByTwisting", {}).get("conformerXyzsForTorsionScans", [])
    if not torsion_scan_xyzs:
        torsion_scan_xyzs = _infer_torsion_scan_xyzs(config, conformer_xyzs)
        if torsion_scan_xyzs:
            config.setdefault("runtimeInfo", {}).setdefault("madeByTwisting", {})[
                "conformerXyzsForTorsionScans"
            ] = torsion_scan_xyzs
    selected_for_charges = set(_conformer_names_from_xyzs(charges_xyzs))
    selected_for_scans = set(_conformer_names_from_xyzs(torsion_scan_xyzs))
    pca_df["selectedForCharges"] = pca_df["conformer"].isin(selected_for_charges)
    pca_df["selectedForTorsionScans"] = pca_df["conformer"].isin(selected_for_scans)
    plot_html = _build_pca_plot_html(pca_df, explained_variance, selected_for_charges, selected_for_scans)
    conformer_html_map = _build_conformer_html_map(
        pca_df["conformer"].tolist(),
        images_dir,
        reporter_dir,
    )

    backbone_aliases = config.get("moleculeInfo", {}).get("backboneAliases")
    ramachandran_plot_html = None
    ramachandran_available = False
    if _backbone_aliases_ready(backbone_aliases):
        phi_psi_dihedrals = config["runtimeInfo"]["madeByTwisting"]["allDihedrals"].get("phiPsiDihedrals", {})
        phi_by_ca, psi_by_ca = _classify_phi_psi_dihedrals(phi_psi_dihedrals, backbone_aliases)
        paired_ca_indexes = sorted(set(phi_by_ca) & set(psi_by_ca))
        phi_psi_pairs = [
            (ca_index, phi_by_ca[ca_index], psi_by_ca[ca_index])
            for ca_index in paired_ca_indexes
        ]
        if phi_psi_pairs:
            rama_rows = _measure_ramachandran_angles(conformer_xyzs, phi_psi_pairs, conformer_energies)
            if rama_rows:
                rama_df = pd.DataFrame(rama_rows)
                ramachandran_plot_html = _build_ramachandran_plot_html(rama_df, selected_for_charges, selected_for_scans)
                ramachandran_available = True

    pca_data = {
        "torsions": torsions,
        "conformerAngles": angle_df.to_dict(orient="records"),
        "pcaCoordinates": pca_df.to_dict(orient="records"),
        "explainedVariance": explained_variance,
        "pcaAvailable": pca_available,
        "plotHtml": plot_html,
        "ramachandranPlotHtml": ramachandran_plot_html,
        "ramachandranAvailable": ramachandran_available,
        "conformerHtmlMap": conformer_html_map,
        "selectedConformersForCharges": sorted(selected_for_charges),
        "selectedConformersForTorsionScans": sorted(selected_for_scans),
        "nConformers": len(conformer_xyzs),
        "nTorsions": len(torsions),
    }

    config["runtimeInfo"]["madeByReporting"]["conformerPcaData"] = pca_data
    return pca_data

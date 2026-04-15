import io
import os.path as p
import re
from typing import Optional

import yaml


PERIODIC_TABLE_MASS = {
    "H": 1.0080,
    "He": 4.0026,
    "Li": 6.9400,
    "Be": 9.0122,
    "B": 10.8100,
    "C": 12.0110,
    "N": 14.0070,
    "O": 15.9990,
    "F": 18.9984,
    "Ne": 20.1797,
    "Na": 22.9898,
    "Mg": 24.3050,
    "Al": 26.9815,
    "Si": 28.0850,
    "P": 30.9738,
    "S": 32.0600,
    "Cl": 35.4500,
    "Ar": 39.9480,
    "K": 39.0983,
    "Ca": 40.0780,
    "Sc": 44.9559,
    "Ti": 47.8670,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.9380,
    "Fe": 55.8450,
    "Co": 58.9332,
    "Ni": 58.6934,
    "Cu": 63.5460,
    "Zn": 65.3800,
    "Ga": 69.7230,
    "Ge": 72.6300,
    "As": 74.9216,
    "Se": 78.9710,
    "Br": 79.9040,
    "Kr": 83.7980,
    "Rb": 85.4678,
    "Sr": 87.6200,
    "Y": 88.9058,
    "Zr": 91.2240,
    "Nb": 92.9064,
    "Mo": 95.9500,
    "Tc": 98.0000,
    "Ru": 101.0700,
    "Rh": 102.9055,
    "Pd": 106.4200,
    "Ag": 107.8682,
    "Cd": 112.4140,
    "In": 114.8180,
    "Sn": 118.7100,
    "Sb": 121.7600,
    "Te": 127.6000,
    "I": 126.9045,
    "Xe": 131.2930,
    "Cs": 132.9055,
    "Ba": 137.3270,
    "La": 138.9055,
    "Ce": 140.1160,
    "Pr": 140.9077,
    "Nd": 144.2420,
    "Pm": 145.0000,
    "Sm": 150.3600,
    "Eu": 151.9640,
    "Gd": 157.2500,
    "Tb": 158.9253,
    "Dy": 162.5000,
    "Ho": 164.9303,
    "Er": 167.2590,
    "Tm": 168.9342,
    "Yb": 173.0540,
    "Lu": 174.9668,
    "Hf": 178.4900,
    "Ta": 180.9479,
    "W": 183.8400,
    "Re": 186.2070,
    "Os": 190.2300,
    "Ir": 192.2170,
    "Pt": 195.0840,
    "Au": 196.9665,
    "Hg": 200.5920,
    "Tl": 204.3800,
    "Pb": 207.2000,
    "Bi": 208.9804,
    "Po": 209.0000,
    "At": 210.0000,
    "Rn": 222.0000,
    "Fr": 223.0000,
    "Ra": 226.0000,
    "Ac": 227.0000,
    "Th": 232.0377,
    "Pa": 231.0359,
    "U": 238.0289,
    "Np": 237.0000,
    "Pu": 244.0000,
    "Am": 243.0000,
    "Cm": 247.0000,
    "Bk": 247.0000,
    "Cf": 251.0000,
    "Es": 252.0000,
    "Fm": 257.0000,
    "Md": 258.0000,
    "No": 259.0000,
    "Lr": 262.0000,
    "Rf": 267.0000,
    "Db": 270.0000,
    "Sg": 269.0000,
    "Bh": 270.0000,
    "Hs": 270.0000,
    "Mt": 278.0000,
    "Ds": 281.0000,
    "Rg": 282.0000,
    "Cn": 285.0000,
    "Nh": 286.0000,
    "Fl": 289.0000,
    "Mc": 290.0000,
    "Lv": 293.0000,
    "Ts": 294.0000,
    "Og": 294.0000,
}

VALID_SECTIONS = {"MASS", "NONBON", "NONB", "MOD4      RE"}
SECTION_BREAKS = {"BOND", "ANGLE", "DIHE", "IMPROPER", "CMAP", "END"}


def _extract_float_tokens(tokens: list[str]) -> list[float]:
    values = []
    for token in tokens:
        cleaned = token.strip().strip(",")
        try:
            values.append(float(cleaned))
        except ValueError:
            continue
    return values


def _split_inline_comment(line: str) -> tuple[str, str]:
    if "!" not in line:
        return line, ""
    idx = line.index("!")
    return line[:idx], line[idx + 1 :].strip()


def _normalize_atom_token(atom_type: str) -> str:
    return re.sub(r"[\+\-0-9*]", "", atom_type)


def _infer_element(atom_type: str, mass: Optional[float]) -> str:
    alpha = _normalize_atom_token(atom_type)
    if not alpha:
        return atom_type[0].upper()

    candidates = []
    if len(alpha) >= 2:
        candidates.append(alpha[:2].capitalize())
    candidates.append(alpha[0].upper())

    known_candidates = [c for c in candidates if c in PERIODIC_TABLE_MASS]
    if not known_candidates:
        return candidates[-1]

    # Use mass to resolve ambiguous atom types like "os" (oxygen sp3) vs Os (osmium).
    if mass is not None:
        best_candidate = min(
            known_candidates,
            key=lambda element: abs(PERIODIC_TABLE_MASS[element] - mass),
        )
        best_delta = abs(PERIODIC_TABLE_MASS[best_candidate] - mass)
        if best_delta <= 3.0:
            return best_candidate

    # Without mass, prefer one-letter fallback (safer for GAFF types like "os" -> O).
    # Preserve true chemical symbols such as "Na", "Cl", etc. if typed as proper symbols.
    if len(alpha) >= 2 and alpha[0].isupper() and alpha[1].islower():
        return known_candidates[0]
    return known_candidates[-1]


def _ensure_atom_entry(force_field_dict: dict, atom_type: str, element: str) -> dict:
    if element not in force_field_dict:
        force_field_dict[element] = {}

    if atom_type not in force_field_dict[element]:
        force_field_dict[element][atom_type] = {
            "MASS": None,
            "CHARGE": None,
            "RADIUS": None,
            "WELL_DEPTH": None,
            "COMMENT": "",
            "ORIGIN_COMMENT": "",
            "MASS_ORIGIN": None,
            "CHARGE_ORIGIN": None,
            "RADIUS_ORIGIN": None,
            "WELL_DEPTH_ORIGIN": None,
            "COMMENT_ORIGIN": None,
        }

    return force_field_dict[element][atom_type]


def _find_existing_entry(force_field_dict: dict, atom_type: str) -> tuple[Optional[str], Optional[dict]]:
    for element, atom_types in force_field_dict.items():
        if atom_type in atom_types:
            return element, atom_types[atom_type]
    return None, None


def _infer_charge_from_type(atom_type: str) -> Optional[float]:
    m = re.search(r"([0-9]+)?([+-])", atom_type)
    if not m:
        return None
    magnitude = float(m.group(1)) if m.group(1) else 1.0
    return -magnitude if m.group(2) == "-" else magnitude


def parse_amber_force_field(file_content: str, source_name: str, constant_charge: bool = False) -> dict:
    """
    Parse AMBER .dat/.frcmod content into:
    {ELEMENT: {ATOM_TYPE: {...parameters..., ...origin comments...}}}
    """
    force_field_dict = {}
    lines = io.StringIO(file_content).readlines()

    section = "MASS"
    seen_explicit_section = False

    for raw_line in lines:
        stripped = raw_line.strip()
        if not stripped:
            continue

        upper = stripped.upper()
        if upper in VALID_SECTIONS:
            section = upper
            seen_explicit_section = True
            continue

        if upper in SECTION_BREAKS:
            section = "NONE"
            continue

        body, inline_comment = _split_inline_comment(stripped)
        parts = body.split()
        if not parts:
            continue

        atom_type = parts[0]
        if not atom_type[0].isalnum():
            continue

        # Old AMBER .dat files often start directly with MASS-like entries.
        active_section = section
        if not seen_explicit_section and active_section not in {"NONBON", "NONB", "MOD4      RE"}:
            active_section = "MASS"

        if active_section == "MASS":
            # In legacy .dat files, bond/angle rows can appear without explicit section headers.
            if "-" in atom_type:
                continue

            floats = _extract_float_tokens(parts[1:])
            if not floats:
                continue

            mass = floats[0]
            charge = None
            if constant_charge:
                charge = _infer_charge_from_type(atom_type)
            elif len(floats) >= 2:
                charge = floats[1]

            text_comment = inline_comment if inline_comment else " ".join(parts[3:]).strip()
            element = _infer_element(atom_type, mass)
            entry = _ensure_atom_entry(force_field_dict, atom_type, element)

            origin = f"{source_name}:MASS"
            entry["MASS"] = mass
            entry["MASS_ORIGIN"] = origin
            if charge is not None:
                entry["CHARGE"] = charge
                entry["CHARGE_ORIGIN"] = origin
            if text_comment:
                entry["COMMENT"] = text_comment
                entry["COMMENT_ORIGIN"] = origin

        elif active_section in {"NONBON", "NONB", "MOD4      RE"}:
            floats = _extract_float_tokens(parts[1:])
            if len(floats) < 2:
                continue

            radius, well_depth = floats[0], floats[1]
            existing_element, entry = _find_existing_entry(force_field_dict, atom_type)
            if entry is None:
                inferred_element = _infer_element(atom_type, None)
                entry = _ensure_atom_entry(force_field_dict, atom_type, inferred_element)
            else:
                inferred_element = existing_element

            if entry["MASS"] is not None:
                better_element = _infer_element(atom_type, entry["MASS"])
                if better_element != inferred_element:
                    force_field_dict.setdefault(better_element, {})[atom_type] = force_field_dict[inferred_element].pop(atom_type)
                    if not force_field_dict[inferred_element]:
                        del force_field_dict[inferred_element]
                    entry = force_field_dict[better_element][atom_type]

            origin = f"{source_name}:{active_section}"
            entry["RADIUS"] = radius
            entry["WELL_DEPTH"] = well_depth
            entry["RADIUS_ORIGIN"] = origin
            entry["WELL_DEPTH_ORIGIN"] = origin

            if inline_comment:
                entry["COMMENT"] = inline_comment
                entry["COMMENT_ORIGIN"] = origin

    for element_data in force_field_dict.values():
        for atom_data in element_data.values():
            origins = [
                atom_data.get("MASS_ORIGIN"),
                atom_data.get("CHARGE_ORIGIN"),
                atom_data.get("RADIUS_ORIGIN"),
                atom_data.get("WELL_DEPTH_ORIGIN"),
                atom_data.get("COMMENT_ORIGIN"),
            ]
            unique_origins = []
            for origin in origins:
                if origin and origin not in unique_origins:
                    unique_origins.append(origin)
            atom_data["ORIGIN_COMMENT"] = "; ".join(unique_origins)

    return force_field_dict


def write_dict_to_yaml(data_dict: dict, file_path: str):
    with open(file_path, "w") as yaml_file:
        yaml.dump(data_dict, yaml_file, sort_keys=False, default_flow_style=False, indent=2)


def merge_and_update_dicts(dicts: list[dict]) -> dict:
    mega_dict = {}
    for param_dict in dicts:
        for element, atom_types in param_dict.items():
            mega_dict.setdefault(element, {})
            mega_dict[element].update(atom_types)
    return mega_dict


def load_force_field_file(file_path: str, constant_charge: bool = False) -> dict:
    with open(file_path, "r") as f:
        file_content = f.read()
    return parse_amber_force_field(file_content=file_content, source_name=p.basename(file_path), constant_charge=constant_charge)


def collect_amber19_nonbonded_parameters(leap_parm_dir: str) -> dict:
    """
    Collect non-bonded lookup data from Amber19-relevant force-field sources
    loaded by AmberTools leaprc files (protein ff19SB + gaff2 + tip3p ions).
    """
    source_specs = [
        ("parm19.dat", False),
        ("frcmod.ff19SB", False),
        ("gaff2.dat", False),
        ("frcmod.tip3p", False),
        ("frcmod.ions1lm_126_tip3p", True),
        ("frcmod.ions234lm_1264_tip3p", True),
    ]

    dicts = []
    for file_name, constant_charge in source_specs:
        file_path = p.join(leap_parm_dir, file_name)
        if not p.exists(file_path):
            continue
        dicts.append(load_force_field_file(file_path, constant_charge=constant_charge))

    return merge_and_update_dicts(dicts)


if __name__ == "__main__":
    leap_parm_dir = "/home/esp/.conda/envs/Igor/dat/leap/parm"
    output_yaml = "/home/esp/drFrankenstein/Laboratory/Experiments/Protocol_4_Assembly/gaff2.yaml"

    force_field_dict = collect_amber19_nonbonded_parameters(leap_parm_dir)
    write_dict_to_yaml(force_field_dict, output_yaml)

from typing import Any


PERIODIC_TABLE: dict[str, int] = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
    "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
    "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109,
    "Ds": 110, "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118,
}


def _normalize_element_symbol(raw_value: Any, atom_name: str) -> str:
    if raw_value is None:
        text = ""
    else:
        text = str(raw_value).strip()

    if not text or text.lower() in {"nan", "none"}:
        atom_text = "".join(char for char in str(atom_name).strip() if char.isalpha())
        if not atom_text:
            raise ValueError("Found an atom with neither a valid ELEMENT field nor a parseable ATOM_NAME.")
        text = atom_text[0]

    if len(text) == 1:
        return text.upper()
    return text[0].upper() + text[1:].lower()


def _validate_atom_element(element_symbol: str, atom_name: str, atom_id: Any) -> None:
    if element_symbol not in PERIODIC_TABLE:
        raise ValueError(
            f"Unknown element '{element_symbol}' for atom '{atom_name}' (ATOM_ID={atom_id}). "
            "Please provide valid PDB element symbols."
        )


def calculate_total_electrons(pdb_df, charge: int) -> int:
    neutral_electron_count = 0
    for _, row in pdb_df.iterrows():
        atom_name = row.get("ATOM_NAME", "?")
        atom_id = row.get("ATOM_ID", "?")
        element_symbol = _normalize_element_symbol(row.get("ELEMENT"), atom_name=atom_name)
        _validate_atom_element(element_symbol, atom_name, atom_id)
        neutral_electron_count += PERIODIC_TABLE[element_symbol]

    total_electrons = neutral_electron_count - charge
    if total_electrons <= 0:
        raise ValueError(
            f"Computed electron count is {total_electrons} after applying charge={charge}. "
            "This is not physically valid."
        )
    return total_electrons


def validate_charge_multiplicity_from_pdb(pdb_df, charge: int, multiplicity: int) -> int:
    if multiplicity < 1:
        raise ValueError(f"Multiplicity must be >= 1, got {multiplicity}.")


    total_electrons = calculate_total_electrons(pdb_df, charge)
    if total_electrons % 2 == multiplicity % 2:
        expected = "odd" if total_electrons % 2 == 0 else "even"
        raise ValueError(
            f"Electron/multiplicity parity mismatch: total electrons={total_electrons} "
            f"({'even' if total_electrons % 2 == 0 else 'odd'}) but multiplicity={multiplicity} "
            f"({'odd' if multiplicity % 2 != 0 else 'even'}). For this electron count, multiplicity "
            f"must be {expected}."
        )

    return total_electrons

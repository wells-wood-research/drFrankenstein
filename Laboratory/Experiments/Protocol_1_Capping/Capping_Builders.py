"""
Capping_Builders.py - Clean Capping Group Builders

This module provides clean builders for NME and ACE capping groups using
proper internal coordinate geometry.

Standard peptide bond parameters:
- C-N bond: 1.33 Å (amide)
- N-H bond: 1.01 Å
- C=O bond: 1.23 Å
- C-C bond: 1.52 Å
- C-N-H angle: ~120°
- C-N-C angle: ~120°
- N-C=O angle: ~120°

Author: drFrankenstein
"""

import numpy as np
import pandas as pd
from typing import Tuple
from . import Capping_Geometry as geom
from . import Capping_Assistant

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

# Standard peptide bond parameters
BOND_C_N = 1.33    # Amide C-N bond
BOND_N_H = 1.01    # N-H bond
BOND_C_O = 1.23    # Carbonyl C=O
BOND_C_C = 1.52    # C-C single bond
ANGLE_PEPTIDE = 120.0  # Standard angle for sp2 amide
DIHEDRAL_TRANS = 180.0  # Trans peptide bond


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def build_nme_cap(mol_df: pd.DataFrame,
                  c_terminal_atom: str,
                  nme_template_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build N-Methyl (NME) capping group at C-terminus using internal coordinates.
    
    NME structure: -C(=O)-NH-CH3
    
    Places atoms in order:
    1. N_N: nitrogen bonded to C-terminus carbonyl
    2. H_N: hydrogen on nitrogen
    3. C_N: methyl carbon
    4. H1_N, H2_N, H3_N: methyl hydrogens (via transformation)
    
    Args:
        mol_df (pd.DataFrame): molecule DataFrame
        c_terminal_atom (str): name of C-terminus atom (carbonyl C)
        nme_template_df (pd.DataFrame): NME template structure
    
    Returns:
        pd.DataFrame: NME cap with properly placed atoms
    """
    nme_df = nme_template_df.copy()
    
    # Get C-terminus coords and bonded atoms
    c_coords = geom.get_coords(mol_df, c_terminal_atom)
    bonded_atoms = Capping_Assistant.find_bonded_atoms(mol_df, c_terminal_atom)
    
    # Find CA (or equivalent) and O atoms
    ca_name = None
    o_name = None
    for atom in bonded_atoms:
        if atom.startswith("C"):
            ca_name = atom
        elif atom.startswith("O"):
            o_name = atom
    
    if not ca_name or not o_name:
        raise ValueError(f"Could not find CA and O atoms bonded to {c_terminal_atom}")
    
    ca_coords = geom.get_coords(mol_df, ca_name)
    o_coords = geom.get_coords(mol_df, o_name)
    
    # 1. Place N_N: nitrogen of NME cap
    # Use sp2 geometry - nitrogen is in plane opposite to O and CA
    nn_coords = geom.place_sp2_substituent(
        center_atom=c_coords,
        bonded_atom1=o_coords,
        bonded_atom2=ca_coords,
        bond_length=BOND_C_N,
        in_plane=True
    )
    nme_df = geom.set_coords(nme_df, "N_N", nn_coords)
    
    # 2. Place H_N: hydrogen on nitrogen
    # Use internal coords: trans to C=O (dihedral ~180°)
    hnn1_coords = geom.place_atom_internal_coords(
        atom_a=ca_coords,
        atom_b=c_coords,
        atom_c=nn_coords,
        bond_length=BOND_N_H,
        angle_deg=ANGLE_PEPTIDE,
        dihedral_deg=DIHEDRAL_TRANS
    )
    nme_df = geom.set_coords(nme_df, "H_N", hnn1_coords)
    
    # 3. Place C_N: methyl carbon
    # Opposite to C and H on the nitrogen (sp2-like geometry)
    cn_coords = geom.place_sp2_substituent(
        center_atom=nn_coords,
        bonded_atom1=c_coords,
        bonded_atom2=hnn1_coords,
        bond_length=BOND_C_C,
        in_plane=True
    )
    nme_df = geom.set_coords(nme_df, "C_N", cn_coords)
    
    # 4. Place methyl hydrogens using SVD alignment
    # Now that we have N_N, H_N, C_N properly placed, align the template
    nme_df = Capping_Assistant.transform_whole(
        originalDf=nme_template_df,
        targetDf=nme_df,
        atomNames=["N_N", "H_N", "C_N"]
    )
    
    return nme_df


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def build_ace_cap(mol_df: pd.DataFrame,
                  n_terminal_atom: str,
                  ace_template_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build Acetyl (ACE) capping group at N-terminus using internal coordinates.
    
    ACE structure: CH3-C(=O)-N-
    
    Places atoms in order:
    1. C_C: carbonyl carbon bonded to N-terminus
    2. O_C: carbonyl oxygen
    3. C2_C: methyl carbon
    4. H1_C, H2_C, H3_C: methyl hydrogens (via transformation)
    
    Args:
        mol_df (pd.DataFrame): molecule DataFrame
        n_terminal_atom (str): name of N-terminus atom
        ace_template_df (pd.DataFrame): ACE template structure
    
    Returns:
        pd.DataFrame: ACE cap with properly placed atoms
    """
    ace_df = ace_template_df.copy()
    
    # Get N-terminus coords and bonded atoms
    n_coords = geom.get_coords(mol_df, n_terminal_atom)
    bonded_atoms = Capping_Assistant.find_bonded_atoms(mol_df, n_terminal_atom)
    
    # Find CA (or equivalent) and H atoms
    ca_name = None
    h_name = None
    for atom in bonded_atoms:
        if atom.startswith("C"):
            ca_name = atom
        elif atom.startswith("H"):
            h_name = atom
    
    # Fallback if no H found (some structures might not have it)
    if not h_name:
        h_name = ca_name
    
    if not ca_name:
        raise ValueError(f"Could not find CA atom bonded to {n_terminal_atom}")
    
    ca_coords = geom.get_coords(mol_df, ca_name)
    if h_name and h_name != ca_name:
        h_coords = geom.get_coords(mol_df, h_name)
    else:
        h_coords = ca_coords
    
    # 1. Place C_C: carbonyl carbon of ACE cap
    # Use sp2-like placement opposite to CA and H
    if h_name != ca_name:
        cc1_coords = geom.place_sp2_substituent(
            center_atom=n_coords,
            bonded_atom1=ca_coords,
            bonded_atom2=h_coords,
            bond_length=BOND_C_N,
            in_plane=True
        )
    else:
        # If no H, use simple direction from CA
        direction = n_coords - ca_coords
        cc1_coords = geom.place_atom_by_bond_length(n_coords, direction, BOND_C_N)
    
    ace_df = geom.set_coords(ace_df, "C_C", cc1_coords)
    
    # 2. Place O_C: carbonyl oxygen
    # Trans to CA (dihedral ~180° for trans peptide)
    oc_coords = geom.place_atom_internal_coords(
        atom_a=ca_coords,
        atom_b=n_coords,
        atom_c=cc1_coords,
        bond_length=BOND_C_O,
        angle_deg=ANGLE_PEPTIDE,
        dihedral_deg=DIHEDRAL_TRANS
    )
    ace_df = geom.set_coords(ace_df, "O_C", oc_coords)
    
    # 3. Place C2_C: methyl carbon
    # Opposite to N and O on the carbonyl carbon (sp2 geometry)
    cc2_coords = geom.place_sp2_substituent(
        center_atom=cc1_coords,
        bonded_atom1=n_coords,
        bonded_atom2=oc_coords,
        bond_length=BOND_C_C,
        in_plane=True
    )
    ace_df = geom.set_coords(ace_df, "C2_C", cc2_coords)
    
    # 4. Place methyl hydrogens using SVD alignment
    ace_df = Capping_Assistant.transform_whole(
        originalDf=ace_template_df,
        targetDf=ace_df,
        atomNames=["C_C", "O_C", "C2_C"]
    )
    
    return ace_df


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def validate_geometry(mol_df: pd.DataFrame,
                      cap_df: pd.DataFrame,
                      terminus_type: str,
                      terminal_atom: str) -> dict:
    """
    Validate capping group geometry by checking bond lengths and angles.
    
    Args:
        mol_df (pd.DataFrame): molecule DataFrame
        cap_df (pd.DataFrame): capping group DataFrame
        terminus_type (str): "N" or "C"
        terminal_atom (str): name of terminal atom
    
    Returns:
        dict: validation results with bond lengths and angles
    """
    results = {}
    
    try:
        if terminus_type == "C":
            # NME cap validation
            c_coords = geom.get_coords(mol_df, terminal_atom)
            nn_coords = geom.get_coords(cap_df, "N_N")
            hnn1_coords = geom.get_coords(cap_df, "H_N")
            cn_coords = geom.get_coords(cap_df, "C_N")
            
            results["C-N_N_bond"] = geom.calculate_distance(c_coords, nn_coords)
            results["N_N-H_N_bond"] = geom.calculate_distance(nn_coords, hnn1_coords)
            results["N_N-C_N_bond"] = geom.calculate_distance(nn_coords, cn_coords)
            results["C-N_N-H_N_angle"] = geom.calculate_angle(c_coords, nn_coords, hnn1_coords)
            results["C-N_N-C_N_angle"] = geom.calculate_angle(c_coords, nn_coords, cn_coords)
            
        elif terminus_type == "N":
            # ACE cap validation
            n_coords = geom.get_coords(mol_df, terminal_atom)
            cc1_coords = geom.get_coords(cap_df, "C_C")
            oc_coords = geom.get_coords(cap_df, "O_C")
            cc2_coords = geom.get_coords(cap_df, "C2_C")
            
            results["N-C_C_bond"] = geom.calculate_distance(n_coords, cc1_coords)
            results["C_C-O_C_bond"] = geom.calculate_distance(cc1_coords, oc_coords)
            results["C_C-C2_C_bond"] = geom.calculate_distance(cc1_coords, cc2_coords)
            results["N-C_C-O_C_angle"] = geom.calculate_angle(n_coords, cc1_coords, oc_coords)
            results["N-C_C-C2_C_angle"] = geom.calculate_angle(n_coords, cc1_coords, cc2_coords)
        
        results["status"] = "valid"
        
    except Exception as e:
        results["status"] = "error"
        results["error"] = str(e)
    
    return results

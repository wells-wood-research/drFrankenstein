"""
Capping_Geometry.py - Internal Coordinate Geometry Functions

This module provides clean, robust geometry functions for placing capping groups
using internal coordinates (bond lengths, angles, dihedrals) rather than 
vector averaging.

Author: drFrankenstein
"""

import numpy as np
import pandas as pd
from typing import Tuple

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def normalize_vector(v: np.ndarray) -> np.ndarray:
    """
    Normalize a vector to unit length
    
    Args:
        v (np.ndarray): vector to normalize
    Returns:
        np.ndarray: normalized vector
    """
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        raise ValueError("Cannot normalize zero vector")
    return v / norm


def get_coords(df: pd.DataFrame, atom_name: str) -> np.ndarray:
    """
    Extract XYZ coordinates for a specific atom from a DataFrame
    
    Args:
        df (pd.DataFrame): molecule DataFrame
        atom_name (str): name of atom
    Returns:
        np.ndarray: [x, y, z] coordinates
    """
    coords = df[df["ATOM_NAME"] == atom_name][["X", "Y", "Z"]].values
    if len(coords) == 0:
        raise ValueError(f"Atom {atom_name} not found in DataFrame")
    return coords[0].astype(float)


def set_coords(df: pd.DataFrame, atom_name: str, coords: np.ndarray) -> pd.DataFrame:
    """
    Set XYZ coordinates for a specific atom in a DataFrame
    
    Args:
        df (pd.DataFrame): molecule DataFrame
        atom_name (str): name of atom
        coords (np.ndarray): [x, y, z] coordinates
    Returns:
        pd.DataFrame: updated DataFrame
    """
    df.loc[df["ATOM_NAME"] == atom_name, ["X", "Y", "Z"]] = coords
    return df


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def place_atom_by_bond_length(anchor_coords: np.ndarray,
                               direction_vector: np.ndarray,
                               bond_length: float) -> np.ndarray:
    """
    Place an atom along a direction vector at a specific distance from anchor
    
    Args:
        anchor_coords (np.ndarray): coordinates of anchor atom
        direction_vector (np.ndarray): direction to place new atom (will be normalized)
        bond_length (float): distance from anchor
    Returns:
        np.ndarray: coordinates of new atom
    """
    unit_direction = normalize_vector(direction_vector)
    return anchor_coords + unit_direction * bond_length


def place_atom_internal_coords(atom_a: np.ndarray,
                                atom_b: np.ndarray,
                                atom_c: np.ndarray,
                                bond_length: float,
                                angle_deg: float,
                                dihedral_deg: float) -> np.ndarray:
    """
    Place a new atom D using internal coordinates (Z-matrix style):
    - Bond D-C with length bond_length
    - Angle D-C-B with angle_deg
    - Dihedral D-C-B-A with dihedral_deg
    
    This is the standard method for building molecules with proper geometry.
    
    Args:
        atom_a (np.ndarray): position of atom A (dihedral reference)
        atom_b (np.ndarray): position of atom B (angle reference)
        atom_c (np.ndarray): position of atom C (bonded atom)
        bond_length (float): C-D bond length in Angstroms
        angle_deg (float): D-C-B angle in degrees
        dihedral_deg (float): D-C-B-A dihedral angle in degrees
    
    Returns:
        np.ndarray: position of new atom D
    """
    # Convert angles to radians
    angle = np.radians(angle_deg)
    dihedral = np.radians(dihedral_deg)
    
    # Vector from C to B
    bc = atom_b - atom_c
    bc = normalize_vector(bc)
    
    # Vector from B to A
    ab = atom_a - atom_b
    ab = normalize_vector(ab)
    
    # Calculate normal to plane ABC (perpendicular vector)
    n = np.cross(ab, bc)
    n_norm = np.linalg.norm(n)
    
    if n_norm < 1e-10:
        # Atoms are collinear, choose arbitrary perpendicular
        # Find a vector not parallel to bc
        if abs(bc[0]) < 0.9:
            perp = np.array([1.0, 0.0, 0.0])
        else:
            perp = np.array([0.0, 1.0, 0.0])
        n = np.cross(bc, perp)
    
    n = normalize_vector(n)
    
    # Vector perpendicular to both n and bc (in the plane)
    m = np.cross(n, bc)
    m = normalize_vector(m)
    
    # Build position using spherical coordinates
    # Start along -bc direction, then rotate by angle and dihedral
    d = bond_length * (
        -bc * np.cos(angle) +
        m * np.sin(angle) * np.cos(dihedral) +
        n * np.sin(angle) * np.sin(dihedral)
    )
    
    return atom_c + d


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def place_sp2_substituent(center_atom: np.ndarray,
                          bonded_atom1: np.ndarray,
                          bonded_atom2: np.ndarray,
                          bond_length: float,
                          in_plane: bool = True) -> np.ndarray:
    """
    Place a substituent on an sp2 hybridized center (e.g., carbonyl carbon).
    
    For sp2 geometry, substituents are in-plane (120° angles) or out-of-plane.
    
    Args:
        center_atom (np.ndarray): sp2 hybridized atom (e.g., C in C=O)
        bonded_atom1 (np.ndarray): first bonded atom (e.g., O)
        bonded_atom2 (np.ndarray): second bonded atom (e.g., N)
        bond_length (float): bond length to new substituent
        in_plane (bool): if True, place in plane; if False, place perpendicular
    
    Returns:
        np.ndarray: position of new substituent
    """
    # Vectors from center to bonded atoms
    v1 = normalize_vector(bonded_atom1 - center_atom)
    v2 = normalize_vector(bonded_atom2 - center_atom)
    
    if in_plane:
        # Place in plane, opposite to average of v1 and v2
        # This gives ~120° angles for sp2 geometry
        direction = -(v1 + v2)
        direction = normalize_vector(direction)
    else:
        # Place perpendicular to plane (cross product)
        direction = np.cross(v1, v2)
        direction = normalize_vector(direction)
    
    return center_atom + direction * bond_length


def place_tetrahedral_substituent(center_atom: np.ndarray,
                                   bonded_atom1: np.ndarray,
                                   bonded_atom2: np.ndarray,
                                   bonded_atom3: np.ndarray,
                                   bond_length: float) -> np.ndarray:
    """
    Place a substituent on a tetrahedral center (sp3 hybridized).
    
    Args:
        center_atom (np.ndarray): tetrahedral atom
        bonded_atom1-3 (np.ndarray): three bonded atoms
        bond_length (float): bond length to new substituent
    
    Returns:
        np.ndarray: position of new substituent
    """
    # Calculate vectors from center to bonded atoms
    v1 = normalize_vector(bonded_atom1 - center_atom)
    v2 = normalize_vector(bonded_atom2 - center_atom)
    v3 = normalize_vector(bonded_atom3 - center_atom)
    
    # Direction opposite to average (tetrahedral geometry)
    direction = -(v1 + v2 + v3)
    direction = normalize_vector(direction)
    
    return center_atom + direction * bond_length


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def calculate_distance(coord1: np.ndarray, coord2: np.ndarray) -> float:
    """Calculate distance between two points"""
    return np.linalg.norm(coord2 - coord1)


def calculate_angle(coord1: np.ndarray, 
                    coord2: np.ndarray, 
                    coord3: np.ndarray) -> float:
    """
    Calculate angle at coord2 (angle 1-2-3) in degrees
    
    Args:
        coord1, coord2, coord3: coordinates of three atoms
    Returns:
        float: angle in degrees
    """
    v1 = coord1 - coord2
    v2 = coord3 - coord2
    
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    # Clamp to avoid numerical errors
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    
    return np.degrees(np.arccos(cos_angle))


def calculate_dihedral(coord1: np.ndarray,
                       coord2: np.ndarray,
                       coord3: np.ndarray,
                       coord4: np.ndarray) -> float:
    """
    Calculate dihedral angle 1-2-3-4 in degrees
    
    Args:
        coord1, coord2, coord3, coord4: coordinates of four atoms
    Returns:
        float: dihedral angle in degrees
    """
    b1 = coord2 - coord1
    b2 = coord3 - coord2
    b3 = coord4 - coord3
    
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    
    n1 = normalize_vector(n1)
    n2 = normalize_vector(n2)
    b2 = normalize_vector(b2)
    
    m1 = np.cross(n1, b2)
    
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    
    return np.degrees(np.arctan2(y, x))


def rotate_points_around_axis(points: np.ndarray,
                              axis_point1: np.ndarray,
                              axis_point2: np.ndarray,
                              angle_deg: float) -> np.ndarray:
    """
    Rotate points around an arbitrary axis using Rodrigues' rotation formula.

    Args:
        points (np.ndarray): N x 3 array of XYZ coordinates.
        axis_point1 (np.ndarray): first point on rotation axis.
        axis_point2 (np.ndarray): second point on rotation axis.
        angle_deg (float): rotation angle in degrees.

    Returns:
        np.ndarray: rotated N x 3 coordinates.
    """
    axis = normalize_vector(axis_point2 - axis_point1)
    angle = np.radians(angle_deg)

    points_centered = points - axis_point1
    cross_term = np.cross(axis, points_centered)
    dot_term = np.dot(points_centered, axis)

    rotated = (
        points_centered * np.cos(angle) +
        cross_term * np.sin(angle) +
        np.outer(dot_term, axis) * (1.0 - np.cos(angle))
    )
    return rotated + axis_point1

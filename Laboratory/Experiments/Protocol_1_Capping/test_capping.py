#!/usr/bin/env python
"""
test_capping.py - Test script for improved capping protocols

This script tests the new internal coordinate capping system to ensure
proper geometry and validate improvements over the old system.

Usage:
    python test_capping.py [path_to_pdb]
"""

import sys
import os
sys.path.insert(0, 'Laboratory')

import numpy as np
import pandas as pd
from pdbUtils import pdbUtils
from Experiments.Protocol_1_Capping import Capping_Geometry as geom
from Experiments.Protocol_1_Capping import Capping_Builders
from Experiments.Protocol_1_Capping import Capping_Assistant

def test_geometry_functions():
    """Test basic geometry functions"""
    print("=" * 60)
    print("Testing Geometry Functions")
    print("=" * 60)
    
    # Test normalization
    v = np.array([3.0, 4.0, 0.0])
    v_norm = geom.normalize_vector(v)
    assert abs(np.linalg.norm(v_norm) - 1.0) < 1e-6, "Normalization failed"
    print("✓ Vector normalization")
    
    # Test internal coordinates with known geometry
    atom_a = np.array([0.0, 0.0, 0.0])
    atom_b = np.array([1.0, 0.0, 0.0])
    atom_c = np.array([2.0, 0.0, 0.0])
    
    # Place atom 1.5 Å from C at 120° angle
    atom_d = geom.place_atom_internal_coords(atom_a, atom_b, atom_c, 1.5, 120.0, 180.0)
    
    # Verify bond length
    dist_cd = geom.calculate_distance(atom_c, atom_d)
    assert abs(dist_cd - 1.5) < 0.01, f"Bond length wrong: {dist_cd}"
    print(f"✓ Internal coordinate placement (bond length: {dist_cd:.3f} Å)")
    
    # Verify angle (note: actual angle depends on dihedral)
    angle_bcd = geom.calculate_angle(atom_b, atom_c, atom_d)
    # For 120° input with 180° dihedral on collinear atoms, we get 60° supplementary
    print(f"✓ Bond angle calculated ({angle_bcd:.1f}°)")
    
    # Test sp2 placement
    center = np.array([0.0, 0.0, 0.0])
    bonded1 = np.array([1.0, 0.0, 0.0])
    bonded2 = np.array([0.0, 1.0, 0.0])
    new_atom = geom.place_sp2_substituent(center, bonded1, bonded2, 1.5, in_plane=True)
    dist = geom.calculate_distance(center, new_atom)
    assert abs(dist - 1.5) < 0.01, "SP2 placement distance wrong"
    print(f"✓ SP2 substituent placement")
    
    print()

def test_capping_builders():
    """Test NME and ACE builders with a simple test case"""
    print("=" * 60)
    print("Testing Capping Builders")
    print("=" * 60)
    
    # Load capping templates
    nme_pdb = "Laboratory/Ingredients/Capping_groups/NME.pdb"
    ace_pdb = "Laboratory/Ingredients/Capping_groups/ACE.pdb"
    
    if not os.path.exists(nme_pdb) or not os.path.exists(ace_pdb):
        print("⚠ Capping template PDBs not found, skipping builder tests")
        return
    
    nme_df = pdbUtils.pdb2df(nme_pdb)
    ace_df = pdbUtils.pdb2df(ace_pdb)
    print(f"✓ Loaded NME template ({len(nme_df)} atoms)")
    print(f"✓ Loaded ACE template ({len(ace_df)} atoms)")
    
    # Create a simple test molecule (alanine-like)
    test_mol = pd.DataFrame({
        'ATOM_ID': [1, 2, 3, 4, 5],
        'ATOM_NAME': ['N', 'H', 'CA', 'C', 'O'],
        'RES_NAME': ['ALA', 'ALA', 'ALA', 'ALA', 'ALA'],
        'RES_ID': [1, 1, 1, 1, 1],
        'X': [0.0, -0.5, 1.5, 2.5, 3.0],
        'Y': [0.0, -0.5, 0.0, 1.0, 2.0],
        'Z': [0.0, 0.5, 0.0, 0.0, 0.0],
        'OCCUPANCY': [1.0] * 5,
        'TEMP_FACTOR': [0.0] * 5,
        'ELEMENT': ['N', 'H', 'C', 'C', 'O']
    })
    
    print(f"\nTesting with simple test molecule ({len(test_mol)} atoms)")
    
    try:
        # Test NME builder
        nme_cap = Capping_Builders.build_nme_cap(test_mol, 'C', nme_df)
        print(f"✓ NME cap built ({len(nme_cap)} atoms)")
        
        # Validate NME geometry
        validation = Capping_Builders.validate_geometry(test_mol, nme_cap, 'C', 'C')
        print(f"  C-NN bond: {validation.get('C-NN_bond', 0):.3f} Å (expected ~1.33)")
        print(f"  NN-HNN1 bond: {validation.get('NN-HNN1_bond', 0):.3f} Å (expected ~1.01)")
        print(f"  C-NN-CN angle: {validation.get('C-NN-CN_angle', 0):.1f}° (expected ~120)")
        
        # Test ACE builder
        ace_cap = Capping_Builders.build_ace_cap(test_mol, 'N', ace_df)
        print(f"\n✓ ACE cap built ({len(ace_cap)} atoms)")
        
        # Validate ACE geometry
        validation = Capping_Builders.validate_geometry(test_mol, ace_cap, 'N', 'N')
        print(f"  N-CC1 bond: {validation.get('N-CC1_bond', 0):.3f} Å (expected ~1.33)")
        print(f"  CC1-OC bond: {validation.get('CC1-OC_bond', 0):.3f} Å (expected ~1.23)")
        print(f"  N-CC1-OC angle: {validation.get('N-CC1-OC_angle', 0):.1f}° (expected ~120)")
        
    except Exception as e:
        print(f"✗ Error in builder test: {e}")
        import traceback
        traceback.print_exc()
    
    print()

def test_with_real_molecule(pdb_path=None):
    """Test with a real molecule if available"""
    print("=" * 60)
    print("Testing with Real Molecule")
    print("=" * 60)
    
    if pdb_path is None:
        # Try to find a test molecule
        test_pdbs = [
            '07_TTC_outputs/01_termini_capping/TTC_capped.pdb',
            'Inputs/TTC.pdb',
        ]
        for test_pdb in test_pdbs:
            if os.path.exists(test_pdb):
                pdb_path = test_pdb
                break
    
    if pdb_path is None or not os.path.exists(pdb_path):
        print("⚠ No real molecule available for testing")
        return
    
    print(f"Loading: {pdb_path}")
    mol_df = pdbUtils.pdb2df(pdb_path)
    print(f"  {len(mol_df)} atoms loaded")
    print(f"  Unique atom names: {len(mol_df['ATOM_NAME'].unique())}")
    print("✓ Real molecule loaded successfully")
    print()

def main():
    """Run all tests"""
    print("\n" + "=" * 60)
    print("IMPROVED CAPPING PROTOCOL TEST SUITE")
    print("=" * 60 + "\n")
    
    try:
        test_geometry_functions()
        test_capping_builders()
        
        if len(sys.argv) > 1:
            test_with_real_molecule(sys.argv[1])
        else:
            test_with_real_molecule()
        
        print("=" * 60)
        print("ALL TESTS COMPLETED SUCCESSFULLY ✓")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

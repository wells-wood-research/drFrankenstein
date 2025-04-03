# BENCHMARKING

We need to decide what methods to use for parameterisation and how to compare them. 
We need a method for:
- Torsion Scanning
- Single-points of torsion scans
- Geometry Optimisation for charge fitting
- Single-points of charge fitting

For Torsion Scanning and Geometry Optimisation methods, we need to compare the geometries that are output:
- RMSD between XYZ coords of scans (hard)
- RMSD between optimised geometries in charge calculations (easy)

For Single-point methods, we need to compare:
- Energies for torsion scans (easy)
- Calculated Charges for charge calculations (easy)

We also need to look at computational cost for each set of methods. Most importantly, we need to assess the scalability of each method set. HF//MP2 may be feasible (if a but slow) for small systems with few rotatable bonds, but will not scale well for larger systems. For these tests, we can look at atom scalability by using molecules with non-rotatable groups (e.g. HIS, PHE, TYR, TRP). We expect the scaling with atoms to depend on the method. 
We expect each method to scale near-linearly with the number of rotatable bonds. We can probe this a series of molecules that have different numbers of rotatable bonds (e.g. GLY, ALA, VAL)

## For TORSIONS
For our torsion scanning, we assessed the accuracy, speed, and scaling of the following methods:

| Opt Method | Single Point Method|
|------------|--------------------|
| XTB2       | MP2                |
| XTB2       | B3LYP-631G(d,p)    |
| XTB2       | revPBE-SVP-D3BJ    |
| HF         | MP2                |
| HF-3c      | B3LYP-631G(d,p)    |
| HF         | MP2                |

Against a small molecule (F-GLY), a medium molecule (F-Phe), and a large molecule (covFAD)
We also investigated the effect of the inclusion of implicit solvation, we expect that each
model will scale as follows with the number of atoms:

| Method            | Gas-Phase Scaling   | Solvent Scaling Impact | Relative Cost Increase |
|-------------------|---------------------|------------------------|------------------------|
| GFN2-xTB          | O(N)–O(N²)          | Stays O(N²)           | 20–50%                |
| revPBE/6-31G      | O(N³)              | Stays O(N³)           | 1.5–3x per cycle      |
| MP2/6-31G         | O(N⁵)              | Stays O(N⁵)           | Minor (~20–40% total) |
| B3LYP/6-31G       | O(N⁴)              | Stays O(N⁴)           | 1.5–3x per cycle      |
| B3LYP/6-311G(d,p) | O(N⁴)              | Stays O(N⁴)           | 1.5–3x per cycle      |

And that the torsion scanning step should scale O(N) with the number of rotatable bonds in the system.


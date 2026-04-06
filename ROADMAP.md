# AMBER COMPATABILITY 
[x] add wildcards to acccount for weird capping groups
    [x] make atomType map using PRMTOP
[x] in the case of re-parameterised Phi/Psi, add extra params to account for ncAA-ncAA interactions

# CHARMM COMPATABILITY
[x] DECL -C and DECL +N
[x] BackBoneAlias["N"] _ C, ""["C"] + N


## CHARGES
[] enforce sum(charges) == totalCharge

# Efficiency
[] problem with df merging in stitching assistant [??]
[] gif creation is too slow (do we turn it off or skip frames?)
[] Stitching loop is mostly waiting for a thread

# debugging
[x] yaml writing after stitching
2. issues with time gantt chart
    [] check how torsion scanning timing is implemented
[] rogue ANTECHAMBER output in cwd
[] wrong labelling of fitting_shuffle PNGs

## CONFIG CONFUSION

[] remove cTermini and nTermini and replace with backboneAliases



## HESSIAN

[x] reorder pipeline IN MAIN
    1. Capping
    2. Wriggle
    3. Charges  [TACK on HESSIAN CALCULATIONS]
    3. Assembly [CLUSTERING BASED ON HESSIAN + CHARGES]
    4. Twist
    6. Stitching 
    7. Creation [CHECK USE OF ATOM TYPES]
    8. Reporter [NEW ATOM TYPE ASSIGNMENT SECTION]

ASSEMBLY PROTOCOL

[x] Find .hess files from charges calculation (might depend on protocol used for charges)
[x] loop through .hess files
    [x] construct hessian matrix from .hess file
    [x] invert to get compliance matrix
    [x] construct adjacancy and angle matrix frpm cappedPdb
    [x] convert adjacency and angle matrix to bond and angle lists
    [x] calculate distance and angle strengths (k values)
    [] extract optimised bond and angle values from cappedPdb

[] average over r0 and k0 values for each bond and angle
[] create a dict with atomIds as keys and:
    [] Element
    [] nBonds
    [] Charge
    [] r0 values
    [] k0 values

STORE A DICT WITH DEFAULT AMBER BACKBONE PARAMS IN IT
AND NON-BONDED AND MASS. 

for each atom type, we need to find the closest gaff/ions atom type (use charge)
and use the NONBONDED params from that atom type.

[] create a script that runs over gaff2.dat, and ions parameters and stores data like this:
```yaml
hc: {MASS: 1.008, NONBONDED: {Radius: 0.0, Well_Depth: 0.0}, CHARGE: 0.25, ELEMENT: H},
hn: {MASS: 1.008, NONBONDED: {Radius: 0.0, Well_Depth: 0.0}, CHARGE: 0.35, ELEMENT: H},
fe2: {MASS: 55.845, NONBONDED: {Radius: 0.0, Well_Depth: 0.0}, CHARGE: 2.0, ELEMENT: Fe},
```
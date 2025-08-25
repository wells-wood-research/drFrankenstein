# AMBER COMPATABILITY 
[x] add wildcards to acccount for weird capping groups
    [x] make atomType map using PRMTOP
[x] in the case of re-parameterised Phi/Psi, add extra params to account for ncAA-ncAA interactions

# CHARMM COMPATABILITY
[] DECL -C and DECL +N
[] BackBoneAlias["N"] _ C, ""["C"] + N


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
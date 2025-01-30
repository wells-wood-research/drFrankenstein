A. **mol2_mm_parameterisation**

Processes mol2 files from the QM scans. They can have GAFF2 atom types, but the CA type should be changed to CX to differentiate it from other sp3 carbons.
Can get frcmod of a reference mol2 file with: 

``parmchk2 -i somemol2filefromQMscans.mol2 -f mol2 -o molecule.frcmod -a Y -p /home/eva/anaconda3/envs/AmberTools22/dat/leap/parm/gaff2.dat``

This frcmod should be updated after a torsion has been fitted and its parameters have been updated.

B. **sire_energy_decomposition**

Calculates MM_total - MM_torsion energy for the torsion being fitted. The resulting array needs to be subtracted from the total QM energy before fitting.

C. **dihedral_fitting**

Will fit all possible cosine functions for a given number of periods on the QM - (MM_total - MM_torsion) data. Gives force constants and phases, which then need to be passed to the frcmod file before processing the mol2 files of the next torsion.

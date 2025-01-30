A. mol2_mm_parameterisation - feed mol2 files from the QM scans. They can have GAFF2 atom types, but the CA type should be changed to CX to differentiate it from other sp3 carbons.
Can get frcmod of a reference mol2 file with: 
parmchk2 -i somemol2filefromQMscans.mol2 -f mol2 -o molecule.frcmod -a Y -p /home/eva/anaconda3/envs/AmberTools22/dat/leap/parm/gaff2.dat

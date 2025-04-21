## CHARMM-STYLE PARAMETERISATION PIPELINE

[x] capping
[x] conformer generation
[x] torison scanning
[] charge calculations
    [x] ORCA SOLVATOR for water interaction energies
    [x] GEOM_OPTS
    [x] SINGLE POINTS
    [x] multiwfn for charge fitting
    [ ] get donor/acceptors
[] parameter fitting
    [] get starting parameter set for molecule
        (CGenFF can be used for this, need a binary license to run)
         
    [] get MM[total] with CHARMM
    [] get mm[torsion] with CHARMM


$K_\chi (1 + \cos(n (\chi - \delta))$

Amplitude * (1 + cos( multiplicty * (Angle - Phase)))


## COMPATIBILITY OPTIONS
[] skipBackboneTorsion
    [] identify backbone torsions? [can we use the graph + N and C termini]

Tried to map capping groups to CHARMM backbone types
PROBLEM Building PSF from resulting RTF and PRM pair
I think the issue is that the CGENFF params for CGENFF--CHARMM params don't exist
SOLUTION:
[] find CGENFF-CGENFF params in par_cgenff.prm
[] copy over to molecule.prm
[] change the CGENFF types to CHARMM types 

## THINGS to loop back into config
[] nShuffles for STITCHING
[] nProcs for GOAT (or use max(16, nCores))
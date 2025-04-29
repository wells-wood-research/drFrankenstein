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
    
[x] parameter fitting
    [x] get starting parameter set for molecule
        (CGenFF can be used for this, need a binary license to run)
         
    [x] get MM[total] with CHARMM
    [x] get mm[torsion] with CHARMM

[] Preservation of backbone parameters
    [] identify BB atoms in moleculeInfo
    [x] get str file, split into RTF and PRM
    [x] copy PRM entries from CGENFF PRM to molecule PRM
    [x] change types of BB + capping groups atoms in RTF
            N -> NH1
            HN -> H
            CA -> CT1
            HA -> HB1
            C -> C
            O -> O
            NN -> NH1
            HNN1 -> H
            CN -> CT1
            HCN1,2,3 -> HB1
            CC1 -> C
            OC -> O
            CC2 -> CT1
            HCC1,2,3 -> HB1

    NOTES 
    So far we have created RTF, PRM and PSF with renamed types. 
    Using Parmed, this was ok. Need to check to see if this will run in openMM
    The RTF is very minimal, and a but useless looking - may be better to ditch it and just use 
    the PRM and PSF together. 


    QUESTIONS:
    [] where do we put this in the pipeline?
        - capping  or parameter fitting?
        - new section for initial params?

## STITCHING
  []  Re-work CHARMM parsing to use PSF files instead of RTF
  []  Use Parmed to update PRM rather than our own parser 



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
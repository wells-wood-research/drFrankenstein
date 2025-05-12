## TIDY UP

[x] make a checkpoint for assembly
[] make a splash for assembly
[] make a splash for drFrankenstein itself
[] edit "what have we created" page to say CHARMM / AMBER
[] build a cleaner
[] make a system for making a PDF
    [] torsion scan energies
    [] charges
    [] parameter fitting
    [] timing

[] add CMAP term for AMBER
[] add optional-ness for backbone preservation for CHARMM
[] delete



## CLEANER

OPTIONS
0: `KEEP ALL`
1: `CAREFUL` - delete only files that are never used and contain no useful information
2: `HARSH` - keep only files that are later used by **drFrankenstein**
3: `BRUTAL` - keep only files mentioned in config 


Files to Keep:

1. Termini Capping
    0: KEEP ALL
    1: KEEP ALL
    2: KEEP ONLY [`MOL_capped.pdb`, `orca_opt.inp`, `orca_opt.out`, `MOL_capped_opt.pdb`]
2. Assembly
    0: KEEP ALL
    1: DELETE `united_capped` DO IN `assembly.py`
    2  KEEP ONLY [`MOL_assembled.*`]
3. GOAT 
    0: KEEP ALL
    1: KEEP [`MOL_conformer_*.xyz`, `GOAT_orca.out`, `GOAT_orca.inp`]
    2: KEEP [`MOL_conformer_*.xyz`, `GOAT_orca.out`, `GOAT_orca.inp`]
4. TWIST
    0: KEEP ALL
    1: REMOVE ORCA JOB files
    2: REMOVE All SP Dirs
5. CHARGES
    0: KEEP ALL
    1: REMOVE ORCA JOB files [QMMM, OPT, SP]
    2: REMOVE SP Dirs
6. PARAMETER FITTING
    0: KEEP ALL
    1: ZIP [`NON-FINAL PARAMS`, `NON-FINAL PNGs`]
    2: KEEP [`FINAL PARAMS`, `FINAL PNGs`, `FINAL GIF`] 
7. CREATION
    0: KEEP ALL
    1: KEEP ALL
    2: KEEP ALL


FOR PARAM FITTING - cope with last PRM / FRCMOD being kept!!!






## LENNARD-JONES PARAMETERS

1. use ORCA SOLVATOR to create a MOLECULE - Helium solvated complex
2. For each ATOM_TYPE:
    For each ATOM_NAME == ATOM_TYPE:
        calculate average distance from atom to nearest He atom
    R_min_RG derived from this.
3. For each ATOM_TYPE:
    For each ATOM_NAME == ATOM_TYPE:
        get nearest He atom
        run a distance scan to get well-depth
    E_min_RG derived from this
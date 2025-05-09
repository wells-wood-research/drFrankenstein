## Parmed implementation

[] use parmed to map backbone defaults to atom types
    [x] CHARMM
    [] AMBER
    [] needs **config** arguments to identify backbone residues 
```yaml
moleculeInfo:
    backboneTypeMap:
        N: ["NX"]
        HN: ["HX"]
        C: ["CX"]
        CA: ["CAX"]
        O: ["OX"]
```

[] use parmed to identify rotatable bonds
    [x] CHARMM
    [] AMBER
    [] re-integrate with twist protocol
        [] so far, just use `uniqueRotatableBonds`
            [] for entries with duplicates, make a **config** option for scanDuplicates (and average) or not
        [] in the future, use `uniqueNonAromaticRotatableBonds`



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
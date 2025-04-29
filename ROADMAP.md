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
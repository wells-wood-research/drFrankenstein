## REPORT

### capping step page ?




## Boltzmann Sampling 
[] Boltazmann-weighted average of energy profiles
[] Boltzmann-weighted sampling for conformers in
    [] torsion scanning
    [] charge fitting

## debug AMBER protocols
[] issues with symmetry and charge fitting - equivalent atoms in different groups??
[x] issues with final creation - wrong config keys?

## making functional params
PROBLEM: BONDS, ANGLES, DIHEDRALS not defined for ncAA -- capping groups
        For now, we are using custom types for capping groups, these need to be the standard UPPERCASE types

SOLUTION:
[x] AMBER: enforce this during the assembly stage [N -> N, C -> C, CA -> CX]
[]  CHARMM: work out whether this is a problem in CHARMM at all

PROBLEM: parameters for interactions between NCAA and adjacent amino acids works, but not when adjacent AA is terminal (XC or CT)

SOLUTION:
[x] AMBER: duplicate all parameters that involve CX, replacing CX with XC and CT
[] CHARMM: work out whether this is a problem in CHARMM at all

## FINAL CREATION
[x] CHARMM: make sure we're getting the right files ie. _51 is not found!


## CAPPING PROTOCOL
[] use backboneAliases to simplify capping protocol
[] remove cTermini and nTermini from config - can use backboneAliases instead

## INSTALLATION
[] test it again on a different server
[] update ORCA to 6.1.0 

## REPORTER
[] check for config change compatability
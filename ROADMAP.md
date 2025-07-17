## REPORT
## Time GANTT
[x] cope with mins/ hours / days for x-axis depending on how long it took

### Torsion Scanning
[x] put atom labels on interactive molecule for torsion atoms

### capping step page

### methods page
[x] create a page
[x] put methods for each section 
[x] gather citations for ORCA
[] hard-code citations for:
    [x] MultiWFN
    [x] RESP
    [x] RESP2
    [x] CGenFF
    [x] Antechamber
    [x] OpenMM
    [>] Various MD methods used in OpenMM single point
        [x] LangevinIntegrator
    [x] Fourier Transform?


[x] put citations on methods page


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
[] CHARMM: make sure we're getting the right files ie. _51 is not found!


## CAPPING PROTOCOL
[] use backboneAliases to simplify capping protocol
[] remove cTermini and nTermini from config - can use backboneAliases instead

## INSTALLATION
[] test it again on a different server
[] update ORCA to 6.1.0 

## REPORTER
[] check for config change compatability
## GENERAL PLAN

[x] PDB input from conformational generation

[] work out what torsions need to be scanned

[x] Generate ORCA input files 

## drCapper.py
Protocol:
1. Load into RDKit [x]
2. Add NMe Ac groups [x]
3. save back to pdb [x]
4. fix pdb col names [x]

## drTwist.py
Protocol: 
1. take pdb file as input [x]          
2. identify torsions that need to be scanned [x]
3. For each torsion that needs to be scanned: 
    While covergance criteria not met:
        a. generate nCpus conformer                                 [x]
        b. perform an opt with xtb                                  [x]
        c. perform a forward scan with xtb                          [x]
        d. perform a backward scan with xtb                         [x]
        e. add results of forwards and backwards scans to dataframe [x]
        f. calculate average                                        [x]
        g. perform convergence check                                [x]
4. plot results [x]
5. return final scan results in a useful format

Things to do:
[x] suppress orca errors
[x] logging
[x] implement while loop
[x] on the fly conformer generation
[] work out what the low-lying scans are (have they secretly exploded?)
[] clean up output directories (more subdirs, delete failed dirs)
[] don't scan amide bonds?
[] suppress "[16:46:20] Molecule does not have explicit Hs. Consider calling AddHs()" warning
[] work out a way of doing (geom-opt) // single-point energy calculations on minima/maxima
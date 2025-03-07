## GENERAL PLAN

01 -> Cap PDB file 
02 -> Generate conformers
03 -> Torsion scanning
05 -> Charge Calculations
06 -> Parameter Fitting 
07 -> Final Creation of Parameter files



## PROTOCOL IMPROVEMENTS

### config checker

### pdb checker
[] enforce that atom names need to be unique in inputPdb
[] RES_NAME cannot be NME or ACE

### drCapper
[x] replace RDKIT with something else 
[] different resIds for capping groups? (this will allow us to not rename atom names for staples etc)
--> pdbUtils style

[x] get N-termini and C-Termini from config["moleculeInfo"]
[x] pre-make PDB files for NME and ACE caps
[x] read inputPdb and acePDB and nmePDB into DataFrames
[x] Delete unwanted atoms bound to N and C termini
[x] align and place cap dfs geometrically - concat
[x] write to cappedPdb


### drTwist

02_torsion_scanning
                | --> conformers [MAKE WITH GOAT]
                | --> torsion_X_Y_Z_A 
                                | --> conformer_1_scans
                                                    | --> optimisation
                                                    | --> forwards_scan
                                                    | --> backwards_scan



[] add to completedScanDirs 
[] look at implicit solvation for XTB scans
[] implement single-point recalculation of All energies
[] implement single-point recalculation of stationary point and mid-points

### drCharge
[] implement RESP2 charge fitting
[] loading bars for multiwfn
[] better regex for multiwfn handling
[] implement user-defined charge groups


## drHybrid
[] in MM_torsion_protocol, antechamber does not like multiple residues

## Miscalaneous tasks
[] work out what the low-lying scans are (have they secretly exploded?)
[] clean up output directories (more subdirs, delete failed dirs)
[] improper terms for Nitrogen centres (1 degree increments)
[] split up pathInfo into better subsections
[] implement SYMMETRY

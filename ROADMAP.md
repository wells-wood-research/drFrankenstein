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
    [] get MM[total] with CHARMM
    [] get mm[torsion] with CHARMM


config thoughts
AMBER options are RESP, RESP2
CHARMM options are SOLVATOR 

We may need to add support foer TIP3P waters instead of ALPB(waters)
This could be done by:
1. SOLVATOR
2. QM/MM opt with TIP3P waters (get params)
3. QM/MM sp with TIP3P waters 
4. MultiWFN charge fitting

we just need to add a different option for SOLVATION in the [chargeFittingProtocol] config field


For the SOLVATOR Method we need the dirs:
04_charge_calculations
    |--> SOLVATOR_calculations
    |--> single_points
    |--> charge_fitting

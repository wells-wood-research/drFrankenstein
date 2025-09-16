# Dr. Frankenstein Configuration File (`config.yaml`) Explained

This document details the parameters within the `config.yaml` file used to control the Dr. Frankenstein molecular parameterization workflow. The configuration is written in YAML format and divided into logical sections.

## Overall Structure
<img src="../Images/Grey_Frank_Basic.png" alt="The Mad Doctor Holds Your Molecule Aloft!" width="200" height="200" style="float: right; margin-right: 15px; margin-bottom: 5px;">

The configuration file is organized into the following main sections:

1.  `pathInfo`: Specifies essential file paths for input, output, and external software.
2.  `moleculeInfo`: Defines the properties and topology of the molecule being parameterized.
3.  `torsionScanInfo`: Configures the torsion scanning procedure.
4.  `chargeFittingInfo`: Configures the charge fitting procedure.
5.  `parameterFittingInfo`: Configures the fitting of torsion parameters to the results of your QM torsion scans.
6.  `miscInfo`: Contains miscellaneous settings for the script execution.
<div style="clear: both;"></div> 

---

## `pathInfo`

This section defines the necessary file system paths for the script to locate input files, external programs, and write output.

*   `inputDir`:
    *   **Description**: The full, absolute path to the directory containing input files (like the initial PDB).
    *   **Type**: `String` (Path)
    *   **Default**: `The current working directory`

*   `outputDir`:
    *   **Description**: The full, absolute path to the directory where all output files will be generated. This directory will be created if it does not exist.
    *   **Type**: `String` (Path)
    *   **Default**: `/$cwd/output` 

*   `multiWfnDir`:
    *   **Description**: The full, absolute path to the installation directory of Multiwfn. The script needs this to find the Multiwfn executable.
    *   **Type**: `String` (Path)
    *   **Default**: `NONE` - absence will throw an error!
*   `orcaExe`:
    *   **Description**: The full, absolute path to the ORCA quantum chemistry software executable.
    *   **Type**: `String` (Path)
    *   **Default**: `NONE` - absence will throw an error!

*   `gaff2Dat`:
    *   **Description**: The full, absolute path to the AMBER `gaff2.dat` parameter file. This is likely used as a reference or starting point, especially if AMBER format output is requested.
    *   **Type**: `String` (Path)
    *   **Default**: automatically finds it using the `$AMBERHOME` environment variable

*   `cgenffExe`:
    *   **Description**: The full, absolute path to the CHARMM's `cgenff` executable. [NOT YET IMPLEMENTED]
    *   **Type**: `String` (Path)
    *   **Default**: `NONE` - without this, the user is propmted to use the `cgenff webserver` instead


---

## `moleculeInfo`

This section provides details about the specific molecule being parameterized.

*   `charge`:
    *   **Description**: The total formal charge of the entire molecule.
    *   **Type**: `Integer`
    *   **Default**: `NONE` - absence will throw an error

*   `multiplicity`:
    *   **Description**: The spin multiplicity of the molecule, calculated as `2S+1` where `S` is the total spin. For typical closed-shell organic molecules, this is `1` (singlet).
    *   **Type**: `Integer` (Usually `1`)
    *   **Default**: `NONE` - absence will throw an error

*   `moleculeName`:
    *   **Description**: A three-letter code used to identify the molecule (residue) in output PDB files and potentially in parameter files. You should have a pdb file with the same name in your input directory.
    *   **Type**: `String` (Typically 3 characters, e.g., `AIB`)
    *   **Default**: `NONE` - absence will throw an error

*   `chargeGroups`:
    *   **Description**: Defines groups of heavy atoms (non-hydrogens) that should have a specific combined integer charge during charge fitting (e.g., RESP). This helps maintain integer charges on chemically relevant functional groups. Any heavy atoms *not* listed in any group will be placed in a "left-overs" group, which will be assigned the remaining charge needed to reach the `moleculeInfo.charge`.
    *   **Type**: `Dictionary`
        *   **Keys**: User-defined names for the charge groups (e.g., `CO`, `CAN`).
        *   **Values**: A dictionary for each group containing:
            *   `atoms`: (`List` of `String`) The names of the heavy atoms belonging to this group.
            *   `charge`: (`Integer`) The target integer charge for this group of atoms.

* `backboneAliases`:
    *   **Description**: A dictionary containing alternative names for backbone atoms (e.g., the main chain nitrogen in an amino acid). Keys in this dictionary must be [`N`, `H` `CA`, `HA` `C`, `O`], corresponding to the main chain nitrogen, its hydrogen, the alpha carbon, its hydrogen, the carbonyl carbon, and the carbonyl oxygen, respectively. In cases where you have a non-canonical amino acid that does not have a standard backbone structure (the chromophore in GFP for example), simply provide aliases for `N` and `C` entries. This will be used for adding capping groups. 
    *   **Type**: `Dictionary`
        *   **Keys**: `N`, `H`, `CA`, `HA`, `C`, `O`
        *   **Values**: `List` of `String` corresponding to the names of the atoms in your molecule.

*Example*
```yaml
## for reparameterising glycine
moleculeInfo:       
  charge: 0                 ## overall charge of 0        
  multiplicity: 1           ## singlet (2S+1)
  moleculeName: GLY         ## three-letter code, there will be a PDB file called GLY.pdb in my input directory
  chargeGroups:             ## charge groups for RESP
    cTerminus:
      atoms: [C, O]         ## carbonyl carbon and carbonyl oxygen set to charge 0
      charge: 0
    nTerminus:
      atoms: [CA, N]        ## alpha carbon and its proton (automatically found!) and main chain nitrogen set to charge 0
      charge: 0
  backboneAliases:          ## backbone aliases
    CA: [CA]        
    HA: [HA]
    C: [C]
    N: [N]
    O: [O]
    H: [HN1]                ## our PDB calls the hydrogen something non-standard!
```



---

## `torsionScanInfo`
<img src="../Images/Grey_Frank_Twist.png" alt="The Mad Doctor Holds Your Molecule Aloft!" width="200" height="200" style="float: right; margin-right: 15px; margin-bottom: 5px;">
This section controls the parameters for the conformational search performed via torsion scanning.

*   `runScansOn`: 
    *   **Description**: A dictionary specifying which types of torsions should be scanned.
    *   **Type**: `Dictionary` - each option is a `Boolean` (True/False)
    *   **Options**:
        * `phiPsi` *(default True)* - setting this to `True` will instruct **drFrankenstein** re-parameterise torsion angles that correspond to the backbone phi (φ) and psi (ψ) dihedral angles. Setting this to `False` will instruct **deFrankenstein** to use the default parameters in for phi and psi angles present in either CHARMM or AMBER force fields (as applicable).
        * `polarProtons` *(default True)* - setting this to `True` will instruct **drFrankenstein** to re-parameterise torsion angles that end in a polar proton. Setting this to `False` skip this step and default parameters will be used.
        * `nonPolarProtons` *(default False)* - setting this to `True` will instruct **drFrankenstein** to re-parameterise torsion angles that end in a non-polar proton. Setting this to `False` skip this step and default parameters will be used.

> NOTE: Setting `polarProtons` and `nonPolarProtons` to `True` will often greatly increase the number of torsions that need to scanned, resulting in a large increase in computational cost.  
*   `nConformers`:
    *   **Description**: The maximum number of low-energy conformers identified during the scan that will be kept for subsequent steps (like single point calculations or charge fitting). Set to `-1` to use all conformers found below a certain energy threshold (threshold not specified here, likely internal).
    *   **Type**: `Integer` (`-1` or positive integer)
    *   **Default**: `-1` (all conformers)

*   `scanMethod`:
    *   **Description**: The quantum mechanics (QM) method used for the geometry optimizations *during* the torsion scan. Typically a computationally cheaper method. Examples: `XTB2`, `PM6`, `revPBE def2-SVP D3BJ`.
    *   **Type**: `String` (ORCA QM method specification)
    *   **Default**: `NONE` - absence will throw an error

*   `scanSolvationMethod`:
    *   **Description**: The implicit solvation model applied during the torsion scan optimizations. Use `Null` or `None` (check script implementation for exact keyword) for gas-phase calculations. Examples: `ALPB(water)`, `CPCM(water)`.
    *   **Type**: `String` (ORCA solvation keyword) or `Null`
    *   **Default**: `Null` - no solvation used!

*   `singlePointMethod`:
    *   **Description**: The QM method used for single point energy calculations on the selected conformers *after* the torsion scan is complete. This is typically a more accurate (and computationally expensive) method than `scanMethod`. Set to `Null` or `None` to skip this step. Examples: `revPBE def2-SVP D3BJ`, `B3LYP/6-31G*`.
    *   **Type**: `String` (ORCA QM method specification) or `Null`
    *   **Default**: `NONE` - absence will throw an error

*   `singlePointSolvationMethod`:
    *   **Description**: The implicit solvation model used for the post-scan single point energy calculations. Use `Null` or `None` for gas-phase calculations.
    *   **Type**: `String` (ORCA solvation keyword) or `Null`
    *   **Default**: `Null` - no solvation used!

**Example**
```yaml
torsionScanInfo:
  runScansOn:
    phiPsi: true                                ## Re-parameterise phi and psi
    polarProtons: false                         ## Skip torsions ending in polar protons
    nonPolarProtons: false                      ## Skip torsions ending in non-polar protons
  nConformers: -1                               ## Use all conformers
  scanMethod: XTB2                              ## Optimise geometries with XTB2 (cheap and fast!)
  scanSolvationMethod: ALPB(water)              ## XTB2 compatible water model
  singlePointMethod: revPBE def2-SVP D3BJ       ## More accurate single point method
  singlePointSolvationMethod: CPCM(water)       ## revPBE compatible water model

```

<div style="clear: both;"></div> 

---

## `chargeFittingInfo`
<img src="../Images/Grey_Frank_Charge.png" alt="The Mad Doctor Holds Your Molecule Aloft!" width="200" height="200" style="float: right; margin-right: 15px; margin-bottom: 5px;">
This section defines the parameters for calculating partial atomic charges.

*   `chargeFittingProtocol`:
    *   **Description**: Specifies the protocol to be used for charge fitting.
    *   **Type**: `String`
    *   **Allowed Values**: `RESP`, `RESP2`, `SOLVATOR` (case-sensitive).

*   `nConformers`:
    *   **Description**: The maximum number of conformers (selected from the torsion scan results or provided separately) to be used in the charge fitting procedure. Using multiple conformers generally leads to more robust charges. Set to `-1` to use all available/selected conformers.
    *   **Type**: `Integer` (`-1` or positive integer)
    *   **Default**: `-1` (all conformers)

*   `nCoresPerCalculation`:
    *   **Description**: The number of CPU cores to allocate for each individual QM calculation performed during the charge fitting process (optimizations and single points).
    *   **Type**: `Integer`
    *   **Default**: `1`

*   `optMethod`:
    *   **Description**: The QM method used for geometry optimization of the selected conformers *before* the final single point calculations used for charge fitting.
    *   **Type**: `String` (ORCA QM method specification)
    *   **Default**: `NONE` - absence will throw an error

*   `optSolvationMethod`:
    *   **Description**: The implicit solvation model applied during the pre-fitting geometry optimizations. Use `Null` or `None` for gas-phase.
    *   **Type**: `String` (ORCA solvation keyword) or `Null`
    *   **Default**: `Null` - no solvation used!

*   `singlePointMethod`:
    *   **Description**: The QM method used for the final single point calculations that generate the electrostatic potential (ESP) data required for charge fitting. This should generally be a reasonably accurate method.
    *   **Type**: `String` (ORCA QM method specification)
    *   **Default**: `NONE` - absence will throw an error

    > NOTE: You cannot use XTB2 for single point calculations! They don't make wavefunction output files, so can't be used for charge fitting.

*   `singlePointSolvationMethod`:
    *   **Description**: The solvation model used for the final single point ESP calculations. This is crucial for obtaining charges appropriate for a specific environment (gas phase or solvent).
    *   **Type**: `String` (ORCA solvation keyword) or `Null`
    *   **Default**: `Null` - no solvation used!

 *   `waterDensity`: 
    *   **Description**: Only used for the `SOLVATOR` charge fitting protocol. This is the number of water molecules to add per nm^2 of molecular surface area. More water molecules added may result in more accurate charges, but at the cost of longer calculation times.
    *   **Type**: `int` (positive number)
    *   **Default**: 10

**Examples**

For AMBER-style charge fitting:
```yaml
chargeFittingInfo:
  chargeFittingProtocol: RESP2              ## will run weighted RESP2 protocol
  nConformers: 10                           ## consider 10 conformers
  nCoresPerCalculation: 1                   ## one core per calculation
  optMethod: XTB2                           ## optimise geometries with cheap XTB2 method
  optSolvationMethod: ALPB(water)           ## XTB2 compatible water model
  singlePointMethod: revPBE def2-SVP D3BJ   ## more accurate single point method
  singlePointSolvationMethod: CPCM(water)   ## revPBE compatible water model
  waterDensity: null                        ## not using SOLVATOR protocol, so leave null
```

For CHARMM-style charge fitting:
```yaml
chargeFittingInfo:
  chargeFittingProtocol: SOLVATOR           ## will run explicit solvent protocol
  nConformers: 10                           ## consider 10 conformers
  nCoresPerCalculation: 1                   ## one core per calculation
  optMethod: XTB2                           ## optimise geometries with cheap XTB2 method
  optSolvationMethod: ALPB(water)           ## XTB2 compatible water model
  singlePointMethod: revPBE def2-SVP D3BJ   ## more accurate single point method
  singlePointSolvationMethod: Null          ## solvent is explicit, so no implicit water needed
  waterDensity: 12                          ## 12 water molecules placed per nm^2 of molecular surface area
```

---

## `parameterFittingInfo`

This section controls the fitting of other bonded and non-bonded parameters, primarily dihedral terms, based on the QM energy profiles obtained from the torsion scan.

*   `forceField`:
    *   **Description**: Specifies the target force field format for the output parameters.
    *   **Type**: `String`
    *   **Allowed Values**: `CHARMM`, `AMBER` (case-sensitive).
    *   **Default**: `NONE` - absence will throw an error

*   `maxCosineFunctions`:
    *   **Description**: The maximum number of cosine terms (multiplicities) to use when fitting the potential energy profile for each scanned dihedral angle.
    *   **Type**: `Integer` (Positive integer, e.g., `4`)
    *   **Default**: `3`

> NOTE:  A higher number allows for more complex energy profiles but increases the risk of over-fitting. 
In our experience, too many terms can lead to unstable fits.
    
*   `maxShuffles`:
    *   **Description**: The maximum number of rounds for the parameter fitting optimization algorithm. If the parameterisation procedure has not converged within this number of iterations, it will terminate.
    *   **Type**: `Integer` (Positive integer, e.g., `10`)
    *   **Default**: `50`

*   `minShuffles`: 
    * **Description**: The minimum number of rounds for the parameter fitting optimization algorithm. Even if the parameterisation procedue has converged, it will not terminate before this number of iterations.
    * **Type**: `Integer` (Positive integer, e.g., `5`)
    * **Default**: `10`

*  `convergenceTolerance`:
    * **Description**: The mean average error (MAE) tolerance (in kcal/mol) for the torsion and total components of the parameter fitting optimization. Used to decide when the parameterisation procedure has converged. 
    * **Type**: `Float` (positive number)
    * **Default**: `Null` - convergence checking disabled

*  `l2DampeningFactor`: 
    * **Description**: The L2 dampening factor for the parameter fitting optimization. This is a regularization term that prevents amplitudes from becoming too large.
    * **Type**: `Float` (positive number)
    * **Default**: `0.1`

*  `sagvolSmoothing`: 
    * **Description**: Whether to apply a SAGVOL smoothing term to the parameter fitting optimization. This can prevent the QM[torsion] term from becoming too choppy, which may lead to unstable fits.
    * **Type**: `Boolean` (`True` or `False`)
    * **Default**: `True`


**Example**
```yaml
parameterFittingInfo:
  forceField: AMBER             ## make AMBER parameters
  maxCosineFunctions: 3         ## use a maximum of 3 cosine terms 
  maxShuffles: 100              ## run for a maximum of 100 iterations
  minShuffles: 20               ## do not converge until at least 20 iterations
  convergenceTolerance: 1.0     ## consider convergence when MAE < 1.0 kcal/mol          
  l2DampingFactor: 0.1          ## apply a L2 dampening factor of 0.1
  sagvolSmoothing: true         ## apply a SAGVOL smoothing term to QM[torsion] term
```
---




## `miscInfo`
<img src="../Images/Grey_Frank_Broom.jpg" alt="The Mad Doctor Holds A Broom!" width="200" height="200" style="float: right; margin-right: 15px; margin-bottom: 5px;">

This section contains miscellaneous settings affecting the overall script execution and resource management.

*   `availableCpus`:
    *   **Description**: The total number of CPU cores available on the machine that the Dr. Frankenstein script can utilize for parallel tasks (like running multiple QM calculations simultaneously). This should not exceed the actual number of cores available.
    *   **Type**: `Integer`
    *   **Default**: Uses `all` available cores
*   `cleanUpLevel`:
    *   **Description**: Controls the deletion of intermediate files generated during the workflow. The higher the value, the more files will be deleted.
    *   **Type**: `int`
    *   **Allowed Values**: `0`, `1`, `2`, `3`
    *   **Default**: `1`

**Example**
```yaml
miscInfo:
  availableCpus: 8      ## use 8 cores
  cleanUpLevel: 2       ## delete intermediate files
```
<div style="clear: both;"></div> 

---
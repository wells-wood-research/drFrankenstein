# Dr. Frankenstein Configuration File (`config.yaml`) Explained

This document details the parameters within the `config.yaml` file used to control the Dr. Frankenstein molecular parameterization workflow. The configuration is written in YAML format and divided into logical sections.

## Overall Structure
<img src="../Images/Grey_Frank_Basic.png" alt="The Mad Doctor Holds Your Molecule Aloft!" width="200" height="200" style="float: right; margin-right: 15px; margin-bottom: 5px;">

The configuration file is organized into the following main sections:

1.  `pathInfo`: Specifies essential file paths for input, output, and external software.
2.  `moleculeInfo`: Defines the properties and topology of the molecule being parameterized.
3.  `torsionScanInfo`: Configures the torsion scanning procedure.
4.  `chargeFittingInfo`: Configures the charge fitting procedure.
5.  `parameterFittingInfo`: Configures the fitting of remaining force field parameters (bonds, angles, torsions).
6.  `miscInfo`: Contains miscellaneous settings for the script execution.
<div style="clear: both;"></div> 

---

## `pathInfo`

This section defines the necessary file system paths for the script to locate input files, external programs, and write output.

*   `inputDir`:
    *   **Description**: The full, absolute path to the directory containing input files (like the initial PDB).
    *   **Type**: `String` (Path)
*   `outputDir`:
    *   **Description**: The full, absolute path to the directory where all output files will be generated. This directory will be created if it does not exist.
    *   **Type**: `String` (Path)
*   `multiWfnDir`:
    *   **Description**: The full, absolute path to the installation directory of Multiwfn. The script needs this to find the Multiwfn executable.
    *   **Type**: `String` (Path)
*   `orcaExe`:
    *   **Description**: The full, absolute path to the ORCA quantum chemistry software executable.
    *   **Type**: `String` (Path)
*   `gaff2Dat`:
    *   **Description**: The full, absolute path to the AMBER `gaff2.dat` parameter file. This is likely used as a reference or starting point, especially if AMBER format output is requested.
    *   **Type**: `String` (Path)

---

## `moleculeInfo`

This section provides details about the specific molecule being parameterized.

*   `charge`:
    *   **Description**: The total formal charge of the entire molecule.
    *   **Type**: `Integer`
*   `multiplicity`:
    *   **Description**: The spin multiplicity of the molecule, calculated as `2S+1` where `S` is the total spin. For typical closed-shell organic molecules, this is `1` (singlet).
    *   **Type**: `Integer` (Usually `1`)
*   `moleculePdb`:
    *   **Description**: The full, absolute path to the input Protein Data Bank (PDB) file containing the initial coordinates of the molecule. The comments suggest it can handle standard amino acid formats (protonated or deprotonated C-terminus).
    *   **Type**: `String` (Path)
*   `moleculeName`:
    *   **Description**: A three-letter code used to identify the molecule (residue) in output PDB files and potentially in parameter files.
    *   **Type**: `String` (Typically 3 characters, e.g., `AIB`)
*   `nTermini`:
    *   **Description**: A list containing the atom name(s) designated as the N-terminus of the molecule (e.g., the main chain nitrogen in an amino acid).
    *   **Type**: `List` of `String`
*   `cTermini`:
    *   **Description**: A list containing the atom name(s) designated as the C-terminus of the molecule (e.g., the main chain carbonyl carbon in an amino acid).
    *   **Type**: `List` of `String`
*   `chargeGroups`:
    *   **Description**: Defines groups of heavy atoms (non-hydrogens) that should have a specific combined integer charge during charge fitting (e.g., RESP). This helps maintain integer charges on chemically relevant functional groups. Any heavy atoms *not* listed in any group will be placed in a "left-overs" group, which will be assigned the remaining charge needed to reach the `moleculeInfo.charge`.
    *   **Type**: `Dictionary`
        *   **Keys**: User-defined names for the charge groups (e.g., `CO`, `CAN`).
        *   **Values**: A dictionary for each group containing:
            *   `atoms`: (`List` of `String`) The names of the heavy atoms belonging to this group.
            *   `charge`: (`Integer`) The target integer charge for this group of atoms.

---

## `torsionScanInfo`
<img src="../Images/Grey_Frank_Twist.png" alt="The Mad Doctor Holds Your Molecule Aloft!" width="200" height="200" style="float: right; margin-right: 15px; margin-bottom: 5px;">
This section controls the parameters for the conformational search performed via torsion scanning.

*   `nConformers`:
    *   **Description**: The maximum number of low-energy conformers identified during the scan that will be kept for subsequent steps (like single point calculations or charge fitting). Set to `-1` to use all conformers found below a certain energy threshold (threshold not specified here, likely internal).
    *   **Type**: `Integer` (`-1` or positive integer)
*   `nScanSteps`:
    *   **Description**: The number of steps (increments) used when scanning each dihedral angle.
    *   **Type**: `Integer`
    *   **Constraint**: The comment `(WARNING, only 37 works at the moment!)` indicates a current limitation. `37` likely corresponds to 10-degree steps (360/10 + 1).
*   `scanMethod`:
    *   **Description**: The quantum mechanics (QM) method used for the geometry optimizations *during* the torsion scan. Typically a computationally cheaper method. Examples: `XTB2`, `PM6`, `GFN2-xTB`.
    *   **Type**: `String` (ORCA QM method specification)
*   `scanSolvationMethod`:
    *   **Description**: The implicit solvation model applied during the torsion scan optimizations. Use `Null` or `None` (check script implementation for exact keyword) for gas-phase calculations. Examples: `ALPB(water)`, `CPCM(water)`.
    *   **Type**: `String` (ORCA solvation keyword) or `Null`
*   `singlePointMethod`:
    *   **Description**: The QM method used for single point energy calculations on the selected conformers *after* the torsion scan is complete. This is typically a more accurate (and computationally expensive) method than `scanMethod`. Set to `Null` or `None` to skip this step. Examples: `revPBE def2-SVP D3BJ`, `B3LYP/6-31G*`.
    *   **Type**: `String` (ORCA QM method specification) or `Null`
*   `singlePointSolvationMethod`:
    *   **Description**: The implicit solvation model used for the post-scan single point energy calculations. Use `Null` or `None` for gas-phase calculations.
    *   **Type**: `String` (ORCA solvation keyword) or `Null`
*   `scanSinglePointsOn`:
    *   **Description**: Specifies which conformers should undergo the single point energy calculation defined above.
    *   **Type**: `String`
    *   **Constraint**: The comment `(WARNING, only works for "all" at the moment!)` indicates that currently, single points must be run on *all* selected conformers if this step is enabled (i.e., if `singlePointMethod` is not `Null`). Allowed value: `all`.
*   `skipPhiPSi`:
    *   **Description**: **(NOT IMPLEMENTED)** Intended to provide an option to exclude the backbone phi (φ) and psi (ψ) dihedral angles from the torsion scan, potentially to preserve secondary structure or improve compatibility with existing force fields.
    *   **Type**: `Boolean` (`true` or `false`)
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
*   `nCoresPerCalculation`:
    *   **Description**: The number of CPU cores to allocate for each individual QM calculation performed during the charge fitting process (optimizations and single points).
    *   **Type**: `Integer`
*   `optMethod`:
    *   **Description**: The QM method used for geometry optimization of the selected conformers *before* the final single point calculations used for charge fitting.
    *   **Type**: `String` (ORCA QM method specification)
*   `optSolvationMethod`:
    *   **Description**: The implicit solvation model applied during the pre-fitting geometry optimizations. Use `Null` or `None` for gas-phase.
    *   **Type**: `String` (ORCA solvation keyword) or `Null`
*   `singlePointMethod`:
    *   **Description**: The QM method used for the final single point calculations that generate the electrostatic potential (ESP) data required for charge fitting. This should generally be a reasonably accurate method.
    *   **Type**: `String` (ORCA QM method specification)
*   `singlePointSolvationMethod`:
    *   **Description**: The solvation model used for the final single point ESP calculations. This is crucial for obtaining charges appropriate for a specific environment (gas phase or solvent).
    *   **Type**: `String` (ORCA solvation keyword) or `Null`
    *   **Special Case**: If using the `SOLVATOR` protocol, this might be set to `TIP3P` or similar to indicate a QM/MM setup with explicit solvent molecules (handled by the SOLVATOR protocol implementation). For standard implicit solvation, use ORCA keywords like `CPCM(water)`.
<div style="clear: both;"></div> 

---

## `parameterFittingInfo`

This section controls the fitting of other bonded and non-bonded parameters, primarily dihedral terms, based on the QM energy profiles obtained from the torsion scan.

*   `forceField`:
    *   **Description**: Specifies the target force field format for the output parameters.
    *   **Type**: `String`
    *   **Allowed Values**: `CHARMM`, `AMBER` (case-sensitive).
*   `maxCosineFunctions`:
    *   **Description**: The maximum number of cosine terms (multiplicities) to use when fitting the potential energy profile for each scanned dihedral angle. A higher number allows for more complex energy profiles but increases the risk of overfitting.
    *   **Type**: `Integer` (Positive integer, e.g., `5`)
*   `nShuffles`:
    *   **Description**: The number of rounds or iterations for the parameter fitting optimization algorithm. More rounds may lead to a better fit but increase computation time.
    *   **Type**: `Integer` (Positive integer, e.g., `10`)

---

## `miscInfo`
<img src="../Images/Grey_Frank_Broom.jpg" alt="The Mad Doctor Holds A Broom!" width="200" height="200" style="float: right; margin-right: 15px; margin-bottom: 5px;">

This section contains miscellaneous settings affecting the overall script execution and resource management.

*   `availableCpus`:
    *   **Description**: The total number of CPU cores available on the machine that the Dr. Frankenstein script can utilize for parallel tasks (like running multiple QM calculations simultaneously). This should not exceed the actual number of cores available.
    *   **Type**: `Integer`
*   `cleanUpLevel`:
    *   **Description**: Controls the deletion of intermediate files generated during the workflow. `basic` likely removes large temporary files while keeping essential logs and results, while `None` keeps everything. `full` and `brutal` would imply more aggressive cleaning (currently not fully implemented).
    *   **Type**: `String`
    *   **Allowed Values**: `None`, `basic`, `full`, `brutal` (case-sensitive).
    *   **Constraint**: The comment `(WARNING: only "None" and "basic" work at the moment!)` indicates a current limitation.
<div style="clear: both;"></div> 

---
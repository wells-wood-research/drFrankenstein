# drFrankenstein Configuration File (`config.yaml`)

This guide matches the current config validator and defaulting logic in `Laboratory/OperatingTools/validate_config.py` and `Laboratory/OperatingTools/set_config_defaults.py`.

The config is split into seven top-level sections:

1. `pathInfo`
2. `moleculeInfo`
3. `torsionScanInfo`
4. `chargeFittingInfo`
5. `parameterFittingInfo`
6. `miscInfo`
7. `conformerGenerationInfo`

The defaults helper also creates a `conformerGenerationInfo` block for conformer-generation settings. This is part of the active config shape used by the workflow and is validated like the other top-level sections.

## Config Helper Ownership

The config logic is split across two helper layers:

* `validate_config.py` checks required structure, types, and allowed values.
* `set_config_defaults.py` injects defaults for missing values and environment-derived paths.

In practice, this means some fields are required by validation, while others are optional but still get sensible defaults when the workflow is run through the defaulting helper first. `conformerGenerationInfo` is one of the sections that both helpers now understand.

## Helper Map

This is the quickest way to remember where each config family is handled:

* `pathInfo`
  * Defaults: `set_config_defaults.py`
  * Validation: `validate_config.py`
  * Used by: all stages that read files, write output, or launch external tools
* `moleculeInfo`
  * Defaults: `set_config_defaults.py`
  * Validation: `validate_config.py`
  * Used by: capping, charging, assembly, torsion generation
* `conformerGenerationInfo`
  * Defaults: `set_config_defaults.py`
  * Validation: `validate_config.py`
  * Used by: `Laboratory/OperatingTools/drOrca.py`
* `torsionScanInfo`
  * Defaults: `set_config_defaults.py`
  * Validation: `validate_config.py`
  * Used by: torsion scanning and stitching
* `chargeFittingInfo`
  * Defaults: `set_config_defaults.py`
  * Validation: `validate_config.py`
  * Used by: charge-fitting protocols
* `parameterFittingInfo`
  * Defaults: `set_config_defaults.py`
  * Validation: `validate_config.py`
  * Used by: torsion fitting and reporting
* `miscInfo`
  * Defaults: `set_config_defaults.py`
  * Validation: `validate_config.py`
  * Used by: conformer selection, logging, assembly choice, runtime controls

## `pathInfo`

Paths to input/output folders and external tools.

### `inputDir`
The directory that contains the starting structure and any other input files.

* Type: `string`
* Required: yes
* Default: none

### `outputDir`
The directory where the workflow writes results.

* Type: `string`
* Required: yes
* Default: none if you validate strictly
* Defaulting behavior: if omitted, `set_config_defaults.py` uses `<cwd>/<moleculeName>_FrankenParams`

### `multiWfnDir`
Path to the Multiwfn installation directory.

* Type: `string`
* Required: yes
* Default: none

### `orcaExe`
Path to the ORCA executable.

* Type: `string`
* Required: yes
* Default: none

### `amberHome`
Path to your AMBER installation.

* Type: `string` or `null`
* Required: no
* Default: `null`
* Defaulting behavior: if omitted, the defaults script tries the `AMBERHOME` environment variable

### `cgenffExe`
Path to the local CGenFF executable.

* Type: `string` or `null`
* Required: no
* Default: `null`
* Meaning of `null`: use the CGenFF server workflow instead of a local executable

### Example
```yaml
pathInfo:
  inputDir: /home/esp/scriptDevelopment/drFrankenstein/Specimens
  outputDir: /home/esp/scriptDevelopment/drFrankenstein/_FrankenParams_AIB
  multiWfnDir: /home/esp/bin/Multiwfn_3.8_dev_bin_Linux_noGUI/
  orcaExe: /home/esp/bin/orca_6_1_0/orca
  amberHome: /home/esp/anaconda3/envs/Igor/
  cgenffExe: null
```

## `moleculeInfo`

Describes the molecule being parameterised.

### `charge`
Total formal charge of the molecule.

* Type: `integer`
* Required: yes
* Default: none

### `multiplicity`
Spin multiplicity, usually `1` for closed-shell molecules.

* Type: `integer`
* Required: yes
* Default: none

### `moleculeName`
Short molecule identifier used in output names and residue labels.

* Type: `string`
* Required: yes
* Default: none

### `chargeGroups`
Optional charge-group definitions used during charge fitting.

* Type: `dictionary` or `null`
* Required: no
* Default: `null`
* Structure: each group must contain:
  * `atoms`: list of atom names
  * `charge`: integer target charge for the group
* Notes:
  * Any atoms not listed in a group are handled as leftovers by the charge-fitting workflow.
  * The validator only checks that each group has the right shape.

### `backboneAliases`
Alternative atom names for the backbone atoms used by the capping and backbone-charge logic.

* Type: `dictionary` or `null`
* Required: no
* Default: `null`
* Expected keys: `N`, `H`, `CA`, `HA`, `C`, `O`
* Values: lists of atom-name strings
* Notes:
  * For terminal capping, the most important entries are `N` and `C`.
  * If `chargeFittingInfo.enforceDefaultBackboneCharges` is `true`, all six keys must be present and each value must be a non-empty list.

### Example
```yaml
moleculeInfo:
  charge: 0
  multiplicity: 1
  moleculeName: AIB
  chargeGroups: null
  backboneAliases:
    N: [N]
    H: [HN]
    CA: [CA]
    HA: [HA]
    C: [C]
    O: [O]
```

## `torsionScanInfo`

Controls relaxed torsion scanning and the conformer pool used to seed it.

### `runScansOn`
Boolean flags that determine which classes of torsions are scanned.

* Type: `dictionary`
* Required: yes
* Default: none
* Keys:
  * `phiPsi`: scan backbone phi/psi torsions
  * `polarProtons`: scan torsions ending in a polar proton
  * `nonPolarProtons`: scan torsions ending in a non-polar proton
  * `amides`: placeholder flag, should remain `false`
  * `nonAromaticRings`: placeholder flag, should remain `false`

### `nConformers`
Maximum number of conformers to pass into torsion scanning.

* Type: `integer`
* Required: yes
* Allowed values: `-1` or any integer `>= 0`
* Meaning of `-1`: use all available conformers

### `scanMethod`
ORCA method used for the torsion-scan geometry optimisations.

* Type: `string`
* Required: yes
* Default: none
* Options: any valid ORCA input method string accepted by your ORCA installation

### `nCoresPerCalculation`
CPU cores assigned to each torsion-scan calculation.

* Type: `integer`
* Required: no in the validator, but defaulted by the helper
* Default: `1`
* Constraint: must be a positive integer, and it cannot exceed `miscInfo.availableCpus`

### `scanSolvationMethod`
Implicit solvent model for the scan optimisations.

* Type: `string` or `null`
* Required: no
* Default: `null`
* Meaning of `null`: gas-phase calculation

### `singlePointMethod`
Optional higher-level ORCA method for post-scan single-point recalculations.

* Type: `string` or `null`
* Required: no
* Default: `null`
* Meaning of `null`: skip the single-point refinement step

### `singlePointSolvationMethod`
Implicit solvent model for post-scan single-point calculations.

* Type: `string` or `null`
* Required: no
* Default: `null`

### Example
```yaml
torsionScanInfo:
  runScansOn:
    phiPsi: true
    polarProtons: false
    nonPolarProtons: false
    amides: false
    nonAromaticRings: false
  nConformers: -1
  scanMethod: XTB2
  nCoresPerCalculation: 1
  scanSolvationMethod: ALPB(water)
  singlePointMethod: null
  singlePointSolvationMethod: null
```

## `chargeFittingInfo`

Controls charge fitting and the QM levels used to generate the electrostatic data.

### `chargeFittingProtocol`
Charge-fitting backend to run.

* Type: `string`
* Required: yes
* Allowed values: `RESP`, `RESP2`, `SOLVATOR`

### `nConformers`
Maximum number of conformers used in charge fitting.

* Type: `integer`
* Required: yes
* Allowed values: `-1` or any integer `>= 0`
* Meaning of `-1`: use all selected conformers

### `nCoresPerCalculation`
CPU cores assigned to each charge-fitting QM calculation.

* Type: `integer`
* Required: yes
* Constraint: must be a positive integer and cannot exceed `miscInfo.availableCpus`

### `optMethod`
ORCA method used for geometry optimisation before the charge-fit calculations.

* Type: `string`
* Required: yes

### `optSolvationMethod`
Implicit solvent model for the optimisation step.

* Type: `string` or `null`
* Required: no
* Default: `null`

### `singlePointMethod`
ORCA method used for the final charge-fitting single-point calculation.

* Type: `string`
* Required: yes

### `singlePointSolvationMethod`
Implicit solvent model used for the final single-point calculation.

* Type: `string` or `null`
* Required: no
* Default: `null`

### `waterDensity`
Used only when `chargeFittingProtocol` is `SOLVATOR`.

* Type: `integer`, `float`, or `null`
* Required: no
* Default: `null`
* Meaning: number of water molecules added per nm^2 of molecular surface area

### `enforceDefaultBackboneCharges`
Forces backbone charge handling to match the default backbone pattern.

* Type: `boolean`
* Required: no
* Default: `false`
* Notes:
  * If `true`, `moleculeInfo.backboneAliases` must contain `N`, `H`, `CA`, `HA`, `C`, and `O`
  * Each of those values must be a non-empty list of strings

### Example
```yaml
chargeFittingInfo:
  chargeFittingProtocol: RESP2
  nConformers: 10
  nCoresPerCalculation: 4
  optMethod: XTB2
  optSolvationMethod: ALPB(water)
  singlePointMethod: revPBE def2-SVP D3BJ
  singlePointSolvationMethod: CPCM(water)
  waterDensity: null
  enforceDefaultBackboneCharges: false
```

## `parameterFittingInfo`

Controls torsion-parameter fitting from the QM scan profiles.

### `forceField`
Target force-field family for the final parameter files.

* Type: `string`
* Required: yes
* Allowed values: `CHARMM`, `AMBER`

### `maxCosineFunctions`
Maximum number of cosine terms used when fitting each torsion.

* Type: `integer`
* Required: yes
* Constraint: must be a positive integer

### `maxShuffles`
Maximum number of fitting rounds before the parameterisation run stops.

* Type: `integer`
* Required: yes
* Constraint: must be a positive integer

### `minShuffles`
Minimum number of fitting rounds before convergence is allowed to stop the run.

* Type: `integer`
* Required: yes
* Constraint: must be a positive integer, and it must be less than or equal to `maxShuffles`

### `l2DampingFactor`
Regularisation term used to damp large Fourier amplitudes.

* Type: `float` or `null`
* Required: no
* Default: `0.1`
* Constraint: if provided as a number, it must be positive

### `sagvolSmoothing`
Whether to apply the SAGVOL smoothing step during torsion fitting.

* Type: `boolean`
* Required: yes
* Default: `true`

### `maeTolTotal`
Legacy/default-only tolerance field retained by the defaults helper.

* Type: `float` or `null`
* Required: no
* Default: `null`
* Notes: present in the defaulting script, but not actively validated in the current schema

### `maeTolTorsion`
Legacy/default-only tolerance field retained by the defaults helper.

* Type: `float` or `null`
* Required: no
* Default: `null`
* Notes: present in the defaulting script, but not actively validated in the current schema

### `converganceTolerance`
Convergence threshold used by the stitching/fitting stage.

* Type: `float` or `null`
* Required: no
* Default: typically `null` unless set by a specimen file
* Notes: the fitting score documentation in `FITTING_SCORE.md` describes how this value is used

### Example
```yaml
parameterFittingInfo:
  forceField: AMBER
  maxCosineFunctions: 3
  maxShuffles: 50
  minShuffles: 10
  l2DampingFactor: 0.1
  sagvolSmoothing: true
  maeTolTotal: null
  maeTolTorsion: null
  converganceTolerance: 0.1
```

## `miscInfo`

General runtime settings, random seeds, and assembly-backend selection.

### `availableCpus`
Total CPU cores available for parallel tasks.

* Type: `integer`
* Required: yes
* Constraint: must be a positive integer

### `cleanUpLevel`
How aggressively the workflow removes intermediate files.

* Type: `integer`
* Required: no
* Default: `1`
* Allowed values: `0`, `1`, `2`, `3`

### `assemblyProtocol`
Parameter-assembly backend to use.

* Type: `string`
* Required: no
* Default: `ANTECHAMBER`
* Allowed values: `ANTECHAMBER`, `CGENFF`, `AGNOSTIC`

### `seed`
Seed used for stochastic conformer selection and fitting order.

* Type: `integer`
* Required: no
* Default: `1818`
* Constraint: must be `0` or greater

### `conformerSelectionMethods`
How conformers are chosen when the workflow needs a subset.

* Type: `string`
* Required: no
* Default: `ENERGY`
* Allowed values: `ENERGY`, `DIVERSE`
* Meaning:
  * `ENERGY`: Boltzmann-weighted selection that prefers low-energy conformers
  * `DIVERSE`: PCA-based clustering that prefers a diverse conformer set

### `debug`
Turns on extra debug behavior in some runtime paths.

* Type: `boolean`
* Required: no
* Default: `false`

### `gpuPlatform`
Preferred OpenMM GPU platform.

* Type: `string` or `null`
* Required: no
* Default: auto-detected by the defaults helper
* Allowed values: `CUDA`, `HIP`, `OpenCL`, `CPU`

### Example
```yaml
miscInfo:
  availableCpus: 8
  cleanUpLevel: 1
  assemblyProtocol: ANTECHAMBER
  seed: 1818
  conformerSelectionMethods: ENERGY
  debug: false
  gpuPlatform: CPU
```

## `conformerGenerationInfo`

Controls the GOAT-based conformer-generation step that seeds later protocols.

### `goatMode`
GOAT search mode.

* Type: `string`
* Required: yes
* Default: `GOAT`
* Common values:
  * `GOAT`
  * `GOAT-ENTROPY`

### `energyCutoff`
Maximum relative energy, in kcal/mol, for keeping a conformer.

* Type: `integer` or `float`
* Required: yes
* Default: `6.0`

### `goatMethod`
Method used by the GOAT conformer search.

* Type: `string`
* Required: yes
* Default: `XTB2`
* Common values:
  * `XTB2`
  * `GFN-FF`

### `conformerSelction`
How a downstream conformer subset is selected from the GOAT output.

* Type: `string`
* Required: yes
* Default: `ENERGY`
* Allowed values: `ENERGY`, `DIVERSE`
* Notes:
  * The key is spelled `conformerSelction` in the current code.
  * `ENERGY` favors low-energy conformers.
  * `DIVERSE` favors conformational diversity.

### Example
```yaml
conformerGenerationInfo:
  goatMode: GOAT-ENTROPY
  energyCutoff: 10.0
  goatMethod: GFN-FF
  conformerSelction: DIVERSE
```

# drFrankenstein Configuration (`config.yaml`)


The `drFrankenstein` configuration file (`config.yaml`) is a YAML file that contains the following sections:

1. `pathInfo`
2. `moleculeInfo`
3. `conformerGenerationInfo`
4. `torsionScanInfo`
5. `chargeFittingInfo`
6. `parameterFittingInfo`
7. `miscInfo`



Several entries in this config require ORCA inputs, these can be found in [ORCA's Manual](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html#hybrid-dft)

Several examples of config/PDB input pairs can be found in the [Specimens Directory](https://github.com/wells-wood-research/drFrankenstein/tree/master/Specimens)

## `pathInfo`

This section lets `drFrankenstein` know about your filesystem.
You will need to tailor this section to your own directory structure.


### `inputDir`
This is the directory containing PDB your input PDB file
- Type: `string`
- Required: yes (after defaults)
- Default: current working directory

### `outputDir`
This is the directory where `drFrankenstein` will write its outputs.
This directory will be created if it does not exist already.
- Type: `string`
- Required: yes (after defaults)
- Default : `<cwd>/<moleculeName>_FrankenParams` when `moleculeName` is valid

### `multiWfnDir`
This is the directory where MultiWFN has been installed
- Type: `string`
- Required: yes
- Validation: must exist and be a directory

### `orcaExe`
This is the path to the ORCA executable within your ORCA install directory
- Type: `string`
- Required: yes
- Validation: must exist and be a file

### `amberHome`
This is the location of your AMBER installation, you can find yours by:
```bash
conda activate Igor
echo $AMBERHOME
```
- Type: `string` or `null`
- Required: no
- Default : value from environment variable `$AMBERHOME` (or left unset if unavailable)




## `moleculeInfo`
This section lets `drFrankenstein` know about your molecule.
You will need to tailor this section to your own molecule.

### `charge`
The net formal charge of the molecule
- Type: `integer`
- Required: yes

### `multiplicity`
The spin multiplicity of the molecule, this will almost always be `1`
- Type: `integer`
- Required: yes

### `moleculeName`
The three-letter name of the molecule.
Your input PDB file will be found at `<inputDir>/<moleculeName>.pdb`
This will be used in the RES_ID column of all PDB files generated (hence the 3-letter name)
- Type: `string`
- Required: yes

### `chargeGroups`

`chargeGroups` are used to constrain fitted RESP charges for adjacent atoms to integer totals.
WARNING: setting your charge groups too small can result in instability in the charge fitting process.
If your parameterisation crashes at the MultiWFN stage, try setting `chargeGroups` to `null` and re-run `drFrankenstein`

Note that only heavy atoms need to be included in `chargeGroups`. To save you some time, `drFrankenstein` will automatically detect adjacent protons and include them in the charge groups.

- Type: `dictionary` or `null`
- Required: no
- Default : `null`
- Structure per group:
  - `atoms`: `list[string]` (required)
  - `charge`: `integer` (required)

Example:
```yaml
chargeGroups:
  backbone:
    atoms: ["N", "CA", "C", "O"]
    charge: 0
  sidechain:
    atoms: ["CB", "CG", "CD", "CE", "CG1", "CD1", "CE1", "CG2", "CD2", "CE2"]
    charge: 0
```

### `backboneAliases`
This section lets `drFrankenstein` which atoms in your molecule act as backbone atoms in a protein chain. 
Normally, the atom names of the backbone atoms are: `N`, `H`, `CA`, `HA`, `C`, `O`. In you molecule, this may not be true. this section helps `drFrankenstein` to map the backbone atoms to the canonical names.
For alpha-amino acids, each each atom in this section should have an entry.
For non-alpha amino acids only the `N` and `C` atoms are needed.
For capping groups, only `N` or `C`  atom mappings are needed.
For ligands, `backboneAliases` should be `null`,. This will skip the capping step. 

Capping groups are added to the atoms identified in `N` and `C` entries. 
Note that each key within `backboneAliases` maps to a `list[string]`. This allows `drFrankenstein` to handle groups that have multiple N and C termini such as stapled peptides.

- Type: `dictionary` or `null`
- Required: no
- Default : `null`
- Extra cross-check when `chargeFittingInfo.enforceDefaultBackboneCharges: true`:
  - Keys `N`, `H`, `CA`, `HA`, `C`, `O` must exist
  - Each value must be a non-empty `list[string]`


## `conformerGenerationInfo`
This section lets `drFrankenstein` know how conformers should be generated.
The options in this section are used to control ORCA's `GOAT` protocol: https://www.faccts.de/docs/orca/6.0/manual/contents/typical/GOAT.html

> NOTE that each option in this section has a default value. You may omit this section entirely and the defaults will be used.



### `goatMode`
This lets you choose between the `GOAT` and `GOAT-ENTROPY` conformer generation methods.
`GOAT` is a global-minimum focused method, while `GOAT-ENTROPY` is a diversity-focused but costlier method.
- Type: `string`
- Required: No
- Allowed values: `GOAT`, `GOAT-ENTROPY`
- Default : `GOAT`

### `energyCutoff`
This sets the energy threshold used to keep/discard conformers. If you observe that all of your conformers look the same, this suggests that your molecule has a particularly "deep" global minimum. In this case, you may want to increase this value.
- Type: `integer` or `float`
- Required: No
- Constraint: must be `>= 0`
- Default : `6.0`


### `goatMethod`
This lets you choose between the `XTB2` and `GFN-FF` conformer generation methods. `XTB2` should produce better conformers, but `GFN-FF` is a lot faster.
- Type: `string`
- Required: No
- Allowed values: `XTB2`, `GFN-FF`
- Default : `XTB2`

### `conformerSelection`
This lets you decide how to select conformers for the charge calculation and torsion scanning protocols.
The `ENERGY` method uses Boltzmann-weighted selection, prioritising conformers with lower energy.
The `DIVERSE` method uses a clustering-based selection strategy, prioritising conformers with higher diversity.
- Type: `string`
- Required: No
- Allowed values: `ENERGY`, `DIVERSE`
- Default : `ENERGY`


## `torsionScanInfo`
This section gives `drFrankenstein` information about the torsion scanning protocol you want to run.

### `runScansOn`
This option is a series of booleans that allows you to skip the scanning/fitting of certain classes of torsion angles. If a torsion is skipped, the by-analogy parameters generated by parmchk2 will be used instead.
While is is possible to parameterise all torsion angles, it can become expensive to do so for large systems. To this end, we provide the option to skip torsion angles with terminal polar [eg. `C-C-N-H`] and non-polar protons [eg. `C-C-C-H`]. 

You may want to use canonical torsion angle parameters to describe your ncAA's Phi and Psi angles. If so, set the `PhiPsi` option to `false`. `drFrankenstein` will automatically assign AMBER19SB's Phi and Psi parameters to your molecule's Phi and Psi angles.


- Type: `dictionary`
- Required: yes
- Validation: each value in this dictionary must be `boolean`
- Defaults helper seeds these keys:
  - `phiPsi` (default `true`)
  - `polarProtons` (default `true`)
  - `nonPolarProtons` (default `false`)


### `nConformers`
This sets the number of conformers to perform torsion scans on. The more conformers you use, the smoother your final QM scan energies should be. However, the more conformers you use, the longer the calculation will take.

- Type: `integer`
- Required: yes
- Constraint: must be `>= -1` (`-1` means use all available)
- Default : `-1`


### `scanMethod`
This is the QM method to be used to perform relaxed torsion scans. This must be a valid main input for ORCA.
Do not include the leading `!` character.
Do not include any solvation keywords (this is handled separately).

- Type: `string`
- Required: yes
- Default: none


### `scanSolvationMethod`
This is the solvation method to be used to perform relaxed torsion scans. This must be a valid solvation method available in ORCA and must be compatible with the `scanMethod` choice. Set to `null` to perform torsion scans in the gas-phase (not recommended).

- Type: `string` or `null`
- Required: no
- Default : `null`



### `singlePointMethod`
This is the QM method to be used to perform optional single-point recalculation of scan windows.
This must be a valid main input for ORCA.
Do not include the leading `!` character.

Set this to `null` to skip single-point recalculation.
- Type: `string` or `null`
- Required: no
- Default : `null`

### `singlePointSolvationMethod`
This is the solvation method to be used to perform optional single-point recalculation of scan windows.
This must be a valid solvation method available in ORCA and must be compatible with the `singlePointMethod` choice.

Set this to `null` to skip single-point recalculation in the gas-phase (not recommended).
- Type: `string` or `null`
- Required: no
- Default : `null`


### `nCoresPerCalculation`
The number of CPUs to use per scan and single-point calculation. If `nConformers` is < your total available CPUs, increase this value, otherwise setting to `1` should be fastest.

- Type: `integer`
- Required: no
- Constraint: must be positive and cannot exceed `miscInfo.availableCpus`
- Default : `1`


## `chargeFittingInfo`
This section gives `drFrankenstein` information about the charge fitting protocol you want to run.


### `chargeFittingProtocol`
This allows you to select either the classic `RESP` charge fitting protocol or the more modern `RESP2` charge fitting protocol.

- Type: `string`
- Required: yes
- Allowed values: `RESP`, `RESP2`


### `nConformers`
This sets the number of conformers to perform charge fitting on. The more more conformers you use, the more representative your charge fitting results will be. However, the more conformers you use, the longer the calculation will take.
- Type: `integer`
- Required: yes
- Constraint: must be `>= -1`
- Default : `-1`


### `nCoresPerCalculation`
The number of CPUs to use per charge fitting calculation. If `nConformers` is < your total available CPUs, increase this value, otherwise setting to `1` should be fastest.
- Type: `integer`
- Required: yes
- Constraint: must be positive
- Default : `1`

### `optMethod`
The QM method to be used to perform the optimization step, prior to the single-point step. This must be a valid main input for ORCA. 
To increase efficiency of the charge fitting protocol, we recommend that you use a lower level of theory for the optimization step and a higher level of theory for the single-point step.
Do not include the leading `!` character.
Do not include any solvation keywords (this is handled separately).
- Type: `string`
- Required: yes

### `optSolvationMethod`
The solvation method to be used to perform the optimization step, prior to the single-point step. This must be a valid solvation method available in ORCA and must be compatible with the `optMethod` choice.
Set to `null` to perform optimization in the gas-phase.
- Type: `string` or `null`
- Required: no
- Default: `null`

### `singlePointMethod`
The QM method to be used to perform the single-point step. This must be a valid main input for ORCA.
Do not include the leading `!` character.
Do not include any solvation keywords (this is handled separately).

- Type: `string`
- Required: yes


### `singlePointSolvationMethod`
The solvation method to be used to perform the single-point step. This must be a valid solvation method available in ORCA and must be compatible with the `singlePointMethod` choice.
Set to `null` to perform single-point in the gas-phase.
If you are using the `RESP2` protocol, you must provide a solvation method here.

- Type: `string` or `null`
- Required: no
- Default: not forced in defaults helper

### `enforceDefaultBackboneCharges`
If set to `true`, `drFrankenstein` will re-set the backbone charges to the default values used in the AMBER19SB forcefield. 
This will only work for alpha-amino acids.

- Type: `boolean`
- Required: no
- Default : `false`
- Cross-section effect: if `true`, strict `moleculeInfo.backboneAliases` checks are enabled.


## `parameterFittingInfo`
This section gives `drFrankenstein` information about the parameter fitting protocol you want to run.
There are a lot of parameters in this section that control how cosine functions are fit to your QM-scan data.
It is hard to tell what values for these parameters are best, we recommend you use the default values to begin with. If your fitting results look bad you can change some values then re-run the fitting by setting `checkpointInfo.torsionFittingComplete` to `false` in your `<outDir>/drFrankenstein.yaml` config file.

### `forceField`
Parameters will be generated in either the AMBER or CHARMM (work in progress) forcefield
- Type: `string`
- Required: yes
- Allowed values: `CHARMM`, `AMBER`

### `maxCosineFunctions`
This sets a hard cap on the number of cosine functions used to fit to each torsion angle scan. An internal protocol will attempt to use as few cosine functions as possible to fit the data. Setting this value too high may cause overfitting, while too low may result in poor descriptions of complex torsional energy landscapes.
- Type: `integer`
- Required: yes
- Constraint: must be positive
- Default : `4`

### `maxShuffles`
This sets a hard cap on the number of "shuffles" (times each torsion is parameterised) used to fit to each torsion angle scan. In most cases, fitting should converge before this limit is reached.
- Type: `integer`
- Required: yes
- Constraint: must be positive
- Default : `50`

### `minShuffles`
This sets the minimum number of "shuffles" (times each torsion is parameterised) before we `drFrankenstein` is allowed to consider a torsion as "converged".
- Type: `integer`
- Required: yes
- Constraint: must be positive and `<= maxShuffles`
- Default : `10`

### `l2DampingFactor`
This sets the L2 dampening factor that gets applied to Amplitudes of cosine functions during the fitting process. Setting this too high will result in parameters with too low energy barriers, while too low can make the fitting process unstable (amplitudes get increasingly large to the point that they are completely unphysical).
- Type: `float` or `null`
- Required: no
- Constraint: if float, must be positive
- Default : `0.1`


### `sagvolSmoothing`
This sets whether or not to perform Savitzky-Golay smoothing during the fitting process.
- Type: `boolean`, `dictionary`, or `null`
- Required: no
- Default : `true`

### `converganceTolerance`
This sets the convergance tolerance for the fitting process. Torsion parameters are considered "converged" if an internal score is less than this value.
We consider a value of `0.1` to be relatively strict, while a value of `0.3` is more is more lenient. 
- Type: `float  `
- Required: yes in validator schema
- Default : `0.1`


## `miscInfo`
This section gives `drFrankenstein` miscellaneous information that cannot be categorised under any other section.

### `availableCpus`
This is a hard-limit on the number of CPUs available to `drFrankenstein`. If left empty, `drFrankenstein` will use all available CPUs.

- Type: `integer`
- Required: yes
- Constraint: must be positive
- Default : system CPU count (`multiprocessing.cpu_count()`)

### `cleanUpLevel`
This sets the level of cleanup performed by `drFrankenstein`. RECOMMENDED VALUE: `1`
- Type: `integer`
- Required: no
- Allowed values: `0..3`
- Default : `1`

### `assemblyProtocol`
This controls how `ATOM_TYPES` are assigned, as well as how bond and angle parameters are generated. 
If set to `ANTECHAMBER`, `drFrankenstein` will use the `ANTECHAMBER` to assign `ATOM_TYPES` and `parmchk2` to generate bond and angle parameters by-analogy to GAFF2.
If set to `AGNOSTIC`, `drFrankenstein` will run frequency calculations during the charge calculation protocol, then obtain bond and angle parameters via hessian matrix inversion. `ATOM_TYPES` are then assigned using a simple K-Means clustering algorithm. Users may wish to use the `AGNOSTIC` protocol in cases where antechamber fails to generate bond and angle parameters.


- Type: `string`
- Required: no
- Allowed values: `ANTECHAMBER`, `AGNOSTIC`
- Default: `ANTECHAMBER`



### `seed`
There are several stochastic components to `drFrankenstein`. This sets the seed used for random number generation, which ensures reproducibility.
- Type: `integer`
- Required: no
- Constraint: must be `>= 0`
- Default : `1818`


### `debug`
This disables paralellisation, which can be useful for debugging.
If you come across any problems, please let us know by submitting an issue on GitHub: `https://github.com/wells-wood-research/drFrankenstein/issues`
- Type: `boolean`
- Required: no
- Default : `false`

### `gpuPlatform`
This sets the GPU platform used by OpenMM during the parameter fitting protocol. Using any GPU platform will speed this process up.
If you don't have access to a GPU (or don't want to use one) set this to `CPU`
- Type: `string` or `null`
- Required: no
- Allowed values: `CUDA`, `HIP`, `OpenCL`, `CPU` (case-insensitive check)
- Default : auto-detected in order `CUDA -> HIP -> OpenCL -> CPU`

## Minimal Example

```yaml
pathInfo:
  inputDir: /home/esp/scriptDevelopment/drFrankenstein/Specimens
  outputDir: /home/esp/scriptDevelopment/drFrankenstein/_FrankenParams_AIB
  multiWfnDir: /home/esp/bin/Multiwfn_3.8_dev_bin_Linux_noGUI/
  orcaExe: /home/esp/bin/orca_6_1_0/orca
  amberHome: /home/esp/anaconda3/envs/Igor/
  cgenffExe: null

moleculeInfo:
  charge: 0
  multiplicity: 1
  moleculeName: AIB
  chargeGroups: null
  backboneAliases: null

conformerGenerationInfo:
  goatMode: GOAT
  energyCutoff: 6.0
  goatMethod: XTB2
  conformerSelection: ENERGY

torsionScanInfo:
  runScansOn:
    phiPsi: true
    polarProtons: true
    nonPolarProtons: false
    amides: false
    nonAromaticRings: false
  nConformers: -1
  scanMethod: XTB2
  nCoresPerCalculation: 1
  scanSolvationMethod: null
  singlePointMethod: null
  singlePointSolvationMethod: null

chargeFittingInfo:
  chargeFittingProtocol: RESP
  nConformers: -1
  nCoresPerCalculation: 1
  optMethod: XTB2
  optSolvationMethod: null
  singlePointMethod: revPBE def2-SVP D3BJ
  singlePointSolvationMethod: null
  waterDensity: null
  enforceDefaultBackboneCharges: false

parameterFittingInfo:
  forceField: AMBER
  maxCosineFunctions: 4
  maxShuffles: 50
  minShuffles: 10
  l2DampingFactor: 0.1
  sagvolSmoothing: true
  maeTolTotal: null
  maeTolTorsion: null
  cosineMinScoreImprovement: 0.0001

miscInfo:
  availableCpus: 8
  cleanUpLevel: 1
  assemblyProtocol: ANTECHAMBER
  seed: 1818
  conformerSelectionMethods: ENERGY
  debug: false
  gpuPlatform: CPU
```

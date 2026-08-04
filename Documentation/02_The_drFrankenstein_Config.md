# drFrankenstein Configuration (`config.yaml`)

The `drFrankenstein` configuration file (`config.yaml`) is a YAML file that contains the following sections:

1. [`pathInfo`](02_The_drFrankenstein_Config_pathInfo.md)
2. [`moleculeInfo`](02_The_drFrankenstein_Config_moleculeInfo.md)
3. [`conformerGenerationInfo`](02_The_drFrankenstein_Config_conformerGenerationInfo.md)
4. [`torsionScanInfo`](02_The_drFrankenstein_Config_torsionScanInfo.md)
5. [`chargeFittingInfo`](02_The_drFrankenstein_Config_chargeFittingInfo.md)
6. [`parameterFittingInfo`](02_The_drFrankenstein_Config_parameterFittingInfo.md)
7. [`miscInfo`](02_The_drFrankenstein_Config_miscInfo.md)

Several entries in this config require ORCA inputs, these can be found in [ORCA's Manual](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html#hybrid-dft)

Several examples of config/PDB input pairs can be found in the [Specimens Directory](https://github.com/wells-wood-research/drFrankenstein/tree/master/Specimens)

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
  debug: false
  gpuPlatform: CPU
```

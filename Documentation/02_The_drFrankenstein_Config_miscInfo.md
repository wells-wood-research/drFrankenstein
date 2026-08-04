# `miscInfo`

This section gives `drFrankenstein` miscellaneous information that cannot be categorised elsewhere.

### `availableCpus`
Hard limit on the CPUs available to `drFrankenstein`. If left empty, all available CPUs are used.
- Type: `integer`
- Required: yes
- Constraint: must be positive
- Default : system CPU count (`multiprocessing.cpu_count()`)

### `cleanUpLevel`
Level of cleanup performed by `drFrankenstein`. Recommended value: `1`.
- Type: `integer`
- Required: no
- Allowed values: `0..3`
- Default : `1`

### `assemblyProtocol`
Controls how `ATOM_TYPES` and bond and angle parameters are generated. `ANTECHAMBER` uses antechamber and parmchk2 by analogy to GAFF2. `AGNOSTIC` obtains bond and angle parameters via Hessian matrix inversion and assigns `ATOM_TYPES` using K-Means clustering.
- Type: `string`
- Required: no
- Allowed values: `ANTECHAMBER`, `AGNOSTIC`
- Default: `ANTECHAMBER`

### `seed`
Seed used for random number generation to ensure reproducibility.
- Type: `integer`
- Required: no
- Constraint: must be `>= 0`
- Default : `1818`

### `debug`
Disables parallelisation, which can be useful for debugging. Problems can be reported at `https://github.com/wells-wood-research/drFrankenstein/issues`.
- Type: `boolean`
- Required: no
- Default : `false`

### `gpuPlatform`
GPU platform used by OpenMM during parameter fitting. Set to `CPU` if no GPU is available.
- Type: `string` or `null`
- Required: no
- Allowed values: `CUDA`, `HIP`, `OpenCL`, `CPU` (case-insensitive check)
- Default : auto-detected in order `CUDA -> HIP -> OpenCL -> CPU`

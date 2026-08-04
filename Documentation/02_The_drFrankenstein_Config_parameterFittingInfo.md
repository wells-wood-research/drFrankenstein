# `parameterFittingInfo`

This section gives `drFrankenstein` information about the parameter fitting protocol. The default values are recommended initially. To re-run fitting after changing values, set `checkpointInfo.torsionFittingComplete` to `false` in `<outDir>/drFrankenstein.yaml`.

### `forceField`
Parameters will be generated in either the AMBER or CHARMM (work in progress) forcefield.
- Type: `string`
- Required: yes
- Allowed values: `CHARMM`, `AMBER`

### `maxCosineFunctions`
Hard cap on the number of cosine functions used for each torsion scan. Too high may cause overfitting; too low may poorly describe complex energy landscapes.
- Type: `integer`
- Required: yes
- Constraint: must be positive
- Default : `4`

### `maxShuffles`
Hard cap on the number of shuffles used to fit each torsion scan.
- Type: `integer`
- Required: yes
- Constraint: must be positive
- Default : `50`

### `minShuffles`
Minimum number of shuffles before `drFrankenstein` considers a torsion converged.
- Type: `integer`
- Required: yes
- Constraint: must be positive and `<= maxShuffles`
- Default : `10`

### `l2DampingFactor`
L2 damping factor applied to cosine amplitudes. Too high gives low energy barriers; too low can make fitting unstable.
- Type: `float` or `null`
- Required: no
- Constraint: if float, must be positive
- Default : `0.1`

### `sagvolSmoothing`
Whether to perform Savitzky-Golay smoothing during fitting.
- Type: `boolean`, `dictionary`, or `null`
- Required: no
- Default : `true`

### `converganceTolerance`
Convergence tolerance for fitting. A value of `0.1` is relatively strict; `0.3` is more lenient.
- Type: `float`
- Required: yes in validator schema
- Default : `0.1`

# `chargeFittingInfo`

This section gives `drFrankenstein` information about the charge fitting protocol you want to run.

### `chargeFittingProtocol`
This allows you to select either the classic `RESP` charge fitting protocol or the more modern `RESP2` charge fitting protocol.
- Type: `string`
- Required: yes
- Allowed values: `RESP`, `RESP2`

### `nConformers`
This sets the number of conformers to perform charge fitting on. More conformers produce more representative results but increase calculation time.
- Type: `integer`
- Required: yes
- Constraint: must be `>= -1`
- Default : `-1`

### `nCoresPerCalculation`
The number of CPUs to use per charge fitting calculation. If `nConformers` is less than your total available CPUs, increase this value; otherwise `1` should be fastest.
- Type: `integer`
- Required: yes
- Constraint: must be positive
- Default : `1`

### `optMethod`
The QM method used for the optimization step before the single-point step. It must be a valid ORCA main input. A lower level of theory is recommended here. Do not include `!` or solvation keywords.
- Type: `string`
- Required: yes

### `optSolvationMethod`
The solvation method used for optimization. It must be compatible with `optMethod`. Set to `null` for gas-phase optimization.
- Type: `string` or `null`
- Required: no
- Default: `null`

### `singlePointMethod`
The QM method used for the single-point step. It must be a valid ORCA main input. Do not include `!` or solvation keywords.
- Type: `string`
- Required: yes

### `singlePointSolvationMethod`
The solvation method used for the single-point step. It must be compatible with `singlePointMethod`. Set to `null` for gas phase. RESP2 requires a solvation method here.
- Type: `string` or `null`
- Required: no
- Default: not forced in defaults helper

### `enforceDefaultBackboneCharges`
If `true`, `drFrankenstein` resets backbone charges to AMBER19SB defaults. This only works for alpha-amino acids.
- Type: `boolean`
- Required: no
- Default : `false`
- Cross-section effect: if `true`, strict `moleculeInfo.backboneAliases` checks are enabled.

# `conformerGenerationInfo`

This section lets `drFrankenstein` know how conformers should be generated. The options control ORCA's `GOAT` protocol: https://www.faccts.de/docs/orca/6.0/manual/contents/typical/GOAT.html

> NOTE that each option in this section has a default value. You may omit this section entirely and the defaults will be used.

### `goatMode`
This lets you choose between the `GOAT` and `GOAT-ENTROPY` conformer generation methods. `GOAT` is a global-minimum focused method, while `GOAT-ENTROPY` is a diversity-focused but costlier method.
- Type: `string`
- Required: No
- Allowed values: `GOAT`, `GOAT-ENTROPY`
- Default : `GOAT`

### `energyCutoff`
This sets the energy threshold used to keep/discard conformers. If all conformers look the same, your molecule may have a particularly deep global minimum; increase this value in that case.
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
This lets you decide how to select conformers for the charge calculation and torsion scanning protocols. `ENERGY` uses Boltzmann-weighted selection, prioritising lower-energy conformers. `DIVERSE` uses clustering, prioritising higher diversity.
- Type: `string`
- Required: No
- Allowed values: `ENERGY`, `DIVERSE`
- Default : `ENERGY`

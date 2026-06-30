# PCY Protocol Test Fixtures

This fixture set is a minimal, fast-testing subset derived from:

- `_FrankenParams_PCY`

It is organized into:

- `inputs/`: representative protocol inputs
- `expected/`: representative expected outputs per protocol stage

The goal is to support protocol-level tests without re-running expensive external jobs
such as ORCA, OpenBabel, tleap, or MultiWFN.

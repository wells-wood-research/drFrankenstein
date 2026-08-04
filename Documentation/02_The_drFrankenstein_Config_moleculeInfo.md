# `moleculeInfo`

This section lets `drFrankenstein` know about your molecule. You will need to tailor this section to your own molecule.

### `charge`
The net formal charge of the molecule
- Type: `integer`
- Required: yes

### `multiplicity`
The spin multiplicity of the molecule, this will almost always be `1`
- Type: `integer`
- Required: yes

### `moleculeName`
The three-letter name of the molecule. Your input PDB file will be found at `<inputDir>/<moleculeName>.pdb`. This will be used in the RES_ID column of all PDB files generated (hence the 3-letter name)
- Type: `string`
- Required: yes

### `chargeGroups`
`chargeGroups` are used to constrain fitted RESP charges for adjacent atoms to integer totals. WARNING: setting your charge groups too small can result in instability in the charge fitting process. If your parameterisation crashes at the MultiWFN stage, try setting `chargeGroups` to `null` and re-run `drFrankenstein`.

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
This section lets `drFrankenstein` know which atoms in your molecule act as backbone atoms in a protein chain. Normally, the atom names of the backbone atoms are: `N`, `H`, `CA`, `HA`, `C`, `O`. In your molecule, this may not be true; this section helps `drFrankenstein` map the backbone atoms to the canonical names.

For alpha-amino acids, each atom in this section should have an entry. For non-alpha amino acids only the `N` and `C` atoms are needed. For capping groups, only `N` or `C` atom mappings are needed. For ligands, `backboneAliases` should be `null`; this will skip the capping step.

Capping groups are added to the atoms identified in `N` and `C` entries. Each key within `backboneAliases` maps to a `list[string]`, allowing `drFrankenstein` to handle groups that have multiple N and C termini such as stapled peptides.
- Type: `dictionary` or `null`
- Required: no
- Default : `null`
- Extra cross-check when `chargeFittingInfo.enforceDefaultBackboneCharges: true`:
  - Keys `N`, `H`, `CA`, `HA`, `C`, `O` must exist
  - Each value must be a non-empty `list[string]`

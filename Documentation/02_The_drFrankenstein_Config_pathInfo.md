# `pathInfo`

This section lets `drFrankenstein` know about your filesystem. You will need to tailor this section to your own directory structure.

### `inputDir`
This is the directory containing PDB your input PDB file
- Type: `string`
- Required: yes (after defaults)
- Default: current working directory

### `outputDir`
This is the directory where `drFrankenstein` will write its outputs. This directory will be created if it does not exist already.
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

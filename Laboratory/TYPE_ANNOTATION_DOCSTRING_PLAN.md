# Laboratory Function Annotation and Docstring Plan

## Goal

Bring every function in `Laboratory/` to a consistent standard:

- all parameters annotated
- all return types annotated
- docstrings present and concise
- no behavioral changes

## Rules

- Do not change logic while doing documentation cleanup.
- Keep docstrings short, consistent, and descriptive.
- Prefer small, reviewable batches over broad edits.
- Verify each batch with `py_compile` before moving on.

## Progress Log

- [x] Shared utility layer: `drFrankenstein.py`, `OperatingTools/drYaml.py`, `OperatingTools/drSubprocess.py`, `OperatingTools/pdb_checker.py`, `OperatingTools/drOrca.py`
- [ ] `OperatingTools/` remaining modules
  - [x] `select_conformers.py`
  - [x] `electron_checker.py`
  - [x] `validate_config.py` docstring pass
- [ ] `Experiments/Protocol_1_Capping`
- [ ] `Experiments/Protocol_2_Wriggling`
- [ ] `Experiments/Protocol_3_Charging`
- [ ] `Experiments/Protocol_4_Assembly`
- [ ] `Experiments/Protocol_5_Twisting`
- [ ] `Experiments/Protocol_6_Stitching`
- [x] `Experiments/Protocol_6_Stitching` core helpers
- [ ] `Experiments/Protocol_7_Creation`
- [x] `Experiments/Protocol_8_Reporter`

## Suggested Batch Order

1. Finish `OperatingTools/`, starting with the most reused modules:
   - `drSplash.py`
   - `select_conformers.py`
   - `drLogger.py`
   - `electron_checker.py`
   - `validate_config.py`
2. Move through `Protocol_8_Reporter`, since it has a lot of helper functions and shared formatting code.
3. Clean up the protocol modules in numerical order from `1` to `7`.

## Working Notes

- If a function already has a good docstring, only add missing type annotations.
- If a function already has annotations, only add or normalize the docstring.
- If a function is private and self-explanatory, keep the docstring brief.
- When a batch is done, mark it complete here and note any edge cases.
- Keep working through the remaining batches until every function under `Laboratory/` is finished; do not stop early just because a subdirectory is complete.

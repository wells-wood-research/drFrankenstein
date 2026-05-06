# drFrankenstein Repo Familiarization Notes

## What this project is
- **drFrankenstein** is an automated molecular-parameter pipeline targeting **AMBER** and **CHARMM** outputs.
- Main run command:
  ```bash
  conda activate Igor
  python ./Laboratory/drFrankenstein.py --config config.yaml
  ```

## High-level pipeline (from `Laboratory/drFrankenstein.py`)
1. Read + default + validate config.
2. Initialize output dir and logger.
3. Run staged protocols with checkpoint gating:
   - Protocol 1: Capping
   - Protocol 2: Conformer generation (Wriggling)
   - Protocol 3: Charge fitting
   - Protocol 4: Parameter assembly
   - Protocol 5: Torsion scanning (Twisting)
   - Protocol 6: Torsion fitting (Stitching)
   - Protocol 7: Final creation
   - Protocol 8: Reporting
4. Write checkpointed config after each stage.

## Key directories/files to know first
- `Laboratory/drFrankenstein.py` – orchestration entry point.
- `Laboratory/Experiments/` – protocol implementations (1..8).
- `Laboratory/OperatingTools/` – shared infra (YAML IO, validation, logging, subprocess helpers, parsers).
- `Documentation/01_Install_Instructions.md` – external deps (ORCA/XTB/MultiWFN, CGenFF notes).
- `Documentation/02_The_drFrankenstein_Config.md` – config schema explanation.
- `Specimens/*.yaml` – example config bundles.
- `tests/` – unit/regression tests.

## Environment and dependencies
- Python target is **3.10**.
- Conda env is defined in `environment.yaml` (`Igor`), with pip extras from `requirements.txt`.
- Critical external tools are expected: **ORCA**, **XTB (otool_xtb for ORCA)**, **MultiWFN**.

## Config shape anchors (quick memory)
- Core sections: `pathInfo`, `moleculeInfo`, `torsionScanInfo`, `chargeFittingInfo`, `parameterFittingInfo`, `miscInfo`.
- `pathInfo.cgenffExe` may be `Null` (server flow) or local executable path.
- `pathInfo.amberHome` is optional.
- `torsionScanInfo.runScansOn` is a boolean flag map.

## Testing quick start
```bash
python -m pytest tests -q
```

Targeted test runs:
```bash
python -m pytest tests/test_laboratory_core.py -q
python -m pytest tests/test_doctors.py -q
python -m pytest tests/test_logging_system.py -q
```

## Logging note (important)
- Code currently initializes logs under:
  - `<outputDir>/00_logs/`
- Some docs mention `<outputDir>/logs/`; treat **code path** as source of truth.

## Fast orientation workflow for future tasks
1. Read `Laboratory/drFrankenstein.py` for stage order and checkpoints.
2. Open the specific protocol doctor module under `Laboratory/Experiments/Protocol_*`.
3. Check shared helper usage under `OperatingTools` before adding new helpers.
4. Validate config assumptions against `validate_config.py` and `set_config_defaults.py`.
5. Run focused tests first, then broader `tests/` as needed.

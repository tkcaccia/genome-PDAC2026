# Shared Runtime

These files were shared across the production pipeline launchers.

- `_common.sh`: common environment loading, runtime selection, and file validation helpers
- `validate_samplesheet.py`: structural samplesheet validator used before each pipeline launch
- `resources.env.example`: sanitized environment template derived from the production workspace

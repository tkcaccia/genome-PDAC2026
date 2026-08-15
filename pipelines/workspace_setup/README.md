# Workspace Setup

This folder contains the machine/bootstrap code used to prepare the production short-read workspace on the remote Ubuntu system.

- `prepare_pancreatic_shortread_workspace.sh`: main end-to-end workspace bootstrap script
- `remote_bootstrap_seq.sh`: remote bootstrap/install helper
- `setup_vep_remote.sh`: offline VEP cache/container setup
- `fix_newvolume_fstab_remote.sh`: persistent mount fix for `/media/user/New Volume3`

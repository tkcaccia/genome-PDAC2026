# Monitoring and Recovery

This folder contains the code used to monitor the remote workstation and automatically resume incomplete PDAC stages.

Files:

- `pdac_remote_watchdog.sh`: remote one-shot/loop watchdog
- `check_pdac_remote_once.sh`: local scheduler helper that syncs and triggers the remote watchdog
- `resume_pdac_queue_remote.sh`: earlier sequential queue wrapper
- `com.stefano.pdac.pipeline.watch.plist`: local `launchd` job used for 30-minute checks
- `ssh_pdac_remote.expect.example`: sanitized example helper; the credentialed production helper is intentionally excluded

#!/usr/bin/env bash
set -Eeuo pipefail

uuid="$(lsblk -no UUID /dev/sda2 | tr -d '[:space:]')"
printf 'UUID=%s\n' "$uuid"

sudo -k cp -a /etc/fstab "/etc/fstab.codex_backup_$(date +%Y%m%d_%H%M%S)"

if ! sudo grep -q "UUID=$uuid" /etc/fstab; then
  printf '\n# Codex PDAC external SSD\nUUID=%s /media/user/New\\040Volume3 ntfs-3g uid=1000,gid=1000,umask=022,nofail,x-gvfs-show 0 0\n' "$uuid" \
    | sudo tee -a /etc/fstab >/dev/null
fi

sudo tail -n 8 /etc/fstab

# Verified MC campaign archive, October 16, 2026

Authorized by the user to copy campaign records to MC_Data_Backup, verify, and then remove verified SSD data. This supersedes the earlier SSD deletion hold.

Source: /Volumes/Office_SSD/Luscher_XY
Destination: /Volumes/MC_Data_Backup/Luscher_XY/archive-2026-10-16

The finite launchd job holds simulation-study, retention-copy, and storage-admission locks. It copies all SSD campaign files and local experiment records plus code, tests, docs and environment specifications, preserving symlinks without following them. The archived chains link points to the archived raw files. All destination regular files are independently reread with SHA-256 before any source deletion. Every SSD source is rehashed before deletion begins; volume UUIDs, closed launch gate, and file identities are checked. Only inventoried SSD files are removed. The live chains symlink is moved to the HDD first. Local study originals remain.

Monitor status.json, transfer.log, checksums.jsonl, verified.json and events.jsonl. A failed_safe state requires inspection; do not claim completion. RunAtLoad is finite, without automatic retry. caffeinate keeps the machine awake while the transfer runs. No simulation kernels or environment changes.

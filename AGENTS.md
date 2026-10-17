# Luscher XY agent instructions

- Work in this checkout: `~/research/Luscher_XY`. The import package is named
  `supervillain`; this is a separate project from `~/research/supervillain_reloaded`.
- Read `PYTHON_ENVIRONMENT.md` for Python setup and verification. Use
  `.venv/bin/python` from the repository root, not bare `python`, system Python,
  or an unqualified `uv run` that could synchronize the environment.
- The current Python/package snapshot is `.python-version` plus
  `requirements-xy.lock.txt`. Do not upgrade, sync, or replace the running
  environment as a side effect of analysis or setup work.
- For XY code changes, run the relevant existing tests; the environment guide
  gives the focused XY suite. Do not launch production campaigns without authorization.
- Inspect `git status --short` before editing. Preserve unrelated modifications
  and untracked campaign files; stage only explicitly intended files.
- Preserve simulation data, raw histories, checkpoints, provenance, and full
  runner logs. Storage cleanup or compaction requires its own explicit scope.
- Read the relevant campaign README and runner/audit configuration before
  interpreting results. Separate correctness checks, descriptive observations,
  and performance claims; report events, denominators, uncertainty, and scope.

# Luscher XY Python environment

The active XY calculation environment uses uv-managed CPython 3.14.6 in
`.venv/`. The macOS system Python is separate and is not used for this project.
`requirements-xy.lock.txt` records all 21 installed distributions as observed
on 2026-09-07. This captures the current environment, not proof of the exact
versions used by every historical run. The snapshot targets macOS arm64 and
pins package versions; it does not pin wheel hashes or every build dependency.

## Use the existing environment

Run from the repository root. No shell activation is required:

```sh
.venv/bin/python your_script.py
.venv/bin/python -m pytest -q test/test_xy.py test/test_xy_cluster.py test/test_xy_scan.py test/test_xy_3d_campaign.py test/test_xy_4d_campaign.py
uv pip check --python .venv/bin/python
```

The test command covers the existing XY action, cluster, checkpoint/replay,
and campaign-audit tests. It is not a production simulation or the full
upstream test suite. Existing tests use temporary output directories.

## Recreate on a fresh checkout

```sh
uv venv --python 3.14.6 .venv
uv pip sync --python .venv/bin/python requirements-xy.lock.txt
```

Use these creation commands only in a fresh checkout, or when deliberately
rebuilding an idle environment. Do not synchronize or replace `.venv` while
calculations are running. Package downloads require network access. Run scripts
from the repository root so they import this checkout's `supervillain` source;
the current environment does not install `supervillain` as a distribution.

## Why there is no new uv.lock

The inherited `pyproject.toml` mixes the numerical dependencies with a much
larger Sphinx/Jupyter toolchain, while `requirements.txt` contains older,
different documentation constraints. The live XY environment has only the
numerical/test dependencies. Running plain `uv sync` or `uv run` can resolve
that broader project specification and change `.venv`; it is not the workflow
for preserving the current calculation environment. Use the explicit commands
above. Separating and validating the documentation toolchain is future work.

No scientific source, existing dependency declaration, or installed package
was changed to introduce this snapshot. Dependency upgrades should be a
separate change with XY validation and a refreshed environment record.

## Setup verification (2026-09-07)

- `uv pip check`: all 21 installed packages compatible.
- Live environment: the focused five-file XY suite passed 162 tests, with
  one skip because its Numba runtime permitted only one thread (12.04 s).
- Recreated temporary environment: all 21 package versions matched the snapshot;
  the same suite passed all 163 tests, including the parallel-thread check (15.34 s).
- The live environment still matched the snapshot after verification. No
  production campaigns or full upstream suite were run.

Test execution used `PYTHONDONTWRITEBYTECODE=1`, temporary `NUMBA_CACHE_DIR`
and `MPLCONFIGDIR`, and pytest `-p no:cacheprovider` so validation did not
write caches into the checkout. Full test logs are temporarily available at
`/private/tmp/luscher-xy-python-setup-tests.log` and
`/private/tmp/luscher-xy-recreated-environment-tests.log`; these are not
permanent research records. The recreated temporary environment was removed
after validation; the live `.venv` was retained.

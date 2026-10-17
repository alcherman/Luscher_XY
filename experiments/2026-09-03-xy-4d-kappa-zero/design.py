"""Freeze the four-dimensional scout manifest."""
import json
from pathlib import Path

from run import OUT, _provenance, atomic_json, now, points, sha256


def main():
    folder = OUT / "batches"
    folder.mkdir(exist_ok=True)
    target = folder / "scout.json"
    manifest = dict(
        name="scout", created_utc=now(), points=points(),
        interpretation="Four evenly spaced epsilon/pi points on the closed interval [1/2,1].",
        update_cycle="one Metropolis sweep, four epsilon-aware Wolff clusters, one winding proposal",
        model_source_sha256=_provenance()["source_sha256"],
        runner_sha256=sha256(OUT / "run.py"),
    )
    if target.exists():
        existing = json.loads(target.read_text())
        for key in ("points", "model_source_sha256", "runner_sha256"):
            if existing[key] != manifest[key]:
                raise RuntimeError(f"Frozen manifest differs in {key}")
        print(f"Verified existing {target}")
    else:
        atomic_json(target, manifest)
        print(f"Wrote {target}")


if __name__ == "__main__":
    main()

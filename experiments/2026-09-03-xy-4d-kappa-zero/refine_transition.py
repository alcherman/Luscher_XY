"""Freeze the second two-point transition-region refinement manifest."""
import json

from run import (OUT, _provenance, atomic_json, now, sha256,
                 transition_refinement_points)


def main():
    folder = OUT / "batches"
    folder.mkdir(exist_ok=True)
    target = folder / "refine-transition.json"
    manifest = dict(
        name="refine-transition", created_utc=now(),
        points=transition_refinement_points(),
        interpretation="Two epsilon/pi values resolving the bracket [0.85, 0.875].",
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

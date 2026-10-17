"""Freeze the three-point refinement manifest inside the scout bracket."""
import json

from run import OUT, _provenance, atomic_json, now, refinement_points, sha256


def main():
    folder = OUT / "batches"
    folder.mkdir(exist_ok=True)
    target = folder / "refine.json"
    manifest = dict(
        name="refine", created_utc=now(), points=refinement_points(),
        interpretation=("Three evenly spaced interior points divide the original "
                        "epsilon/pi bracket [5/6,1] into four equal cells."),
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

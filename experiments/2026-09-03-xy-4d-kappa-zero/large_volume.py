"""Freeze the L=20 check at the two near-critical bracket endpoints."""
import json

from run import OUT, _provenance, atomic_json, large_volume_points, now, sha256


def main():
    folder = OUT / "batches"
    folder.mkdir(exist_ok=True)
    target = folder / "large-volume.json"
    manifest = dict(
        name="large-volume", created_utc=now(), points=large_volume_points(),
        interpretation="L=20 finite-size-flow check at epsilon/pi=0.85 and 7/8.",
        selection_basis="L=8,12,16 place 0.85 on the ordered side and 7/8 on the disordered side.",
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

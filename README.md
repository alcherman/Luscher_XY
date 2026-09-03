# supervillain

[![DOI](https://zenodo.org/badge/679369801.svg)](https://zenodo.org/badge/latestdoi/679369801)

Supervillain is a python package for studying the Villain model.

This fork also supports **D-dimensional XY models with Lüscher admissibility**,
with Python oracles, serial/parallel Numba Metropolis samplers, an XY-only
admissibility-aware Wolff cluster update, vortex and winding
diagnostics, and resumable parameter scans. See the [XY simulation guide](docs/xy.md)
for the action convention, examples, and validation limits.

# Installation + Development

Navigate to the cloned repo and try

```
pip install .  # for development use pip install -e . 
./test/end-to-end.py
```

If pip installation succeeds so too should the example script, which by default samples the (φ, n) formulation of the model on a small lattice.

supervillain has documentation built with [sphinx](https://www.sphinx-doc.org/en/master/).
To build the documentation once you 

```
sphinx-build . _build
```

and then can open `_build/index.html` in a browser.
If you are developing you can replace `sphinx-build` with [`sphinx-autobuild`](https://pypi.org/project/sphinx-autobuild/) to get live updates to the documentation as you go.

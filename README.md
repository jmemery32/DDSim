# DDSim — Damage and Durability Simulator (Level I)

DDSim Level I is a hierarchical, probabilistic fatigue-life predictor: given a
linear-elastic finite-element stress field, it seeds an elliptical flaw at
each candidate node (or a Monte Carlo set of flaws, from an analytical
distribution or a microstructural particle-cracking filter), grows the crack
under NASGRO/Willenborg fatigue crack growth with constant- or
variable-amplitude loading, and reports the predicted life. It was developed
for the dissertation

> J.M. Emery, J.D. Hochhalter, P.A. Wawrzynek, G. Heber, A.R. Ingraffea,
> "DDSim: A hierarchical, probabilistic, multiscale damage and durability
> simulation system – Part I: Methodology and Level I," *Engineering Fracture
> Mechanics* 76 (2009) 1500–1530.

This repository is a from-scratch Python 3 port of the original 2007
Windows/Python 2.4 codebase (which relied on several compiled C++ Python
extensions). See `docs/PORTING_NOTES.md` for the full history of the port —
what changed, what was verified against the original, and every known
behavioral difference.

## Status

**Level I is ported, numba-accelerated, and tested** (178 tests: unit tests
for every module, golden-output regression against the *actual captured
stdout of the original 2006 program*, and a validation suite against the
published dissertation's real 140k-element coupon model). Levels II and III
were incomplete prototypes in the original work and are kept, untouched and
unmaintained, in `legacy/level2_level3_prototype/` — not part of the
supported package.

Not yet done: an Exodus II mesh/stress reader (the current reader is the
original ASCII "RDB" format: `.con`/`.nod`/`.sig`/`.smp`/`.edg`), and
multiprocessing to replace the original Windows-cluster/MPI parallel workflow
(`legacy/windows_cluster_scripts/`).

## Install

```bash
git clone git@github.com:jmemery32/DDSim.git
cd DDSim
make env              # conda env create/update + editable pip install of ddsim
conda activate ddsim
```

`make env` needs conda (or miniconda/mamba) on `PATH` and runs
`conda env update -f environment.yml -n ddsim`, which creates the `ddsim`
conda environment (Python, numpy, scipy, numba, netCDF4, matplotlib, pytest)
and then `pip install`s this repository into it in editable mode.

Already have a Python environment with numpy/scipy/numba? Skip conda and just run:

```bash
pip install -e ".[dev]"
```

## Quickstart

```bash
make example
```

runs the bundled `example1` case (a 20×20×20 cube under uniform stress) for
two nodes. Equivalently, directly:

```bash
ddsim -base example1 -conpath examples/example1/ -parpath examples/example1/ \
      -v -doid_list 10,0 -scale 100
```

`ddsim -help` lists the command-line switches (input file selection, Monte
Carlo vs. deterministic, constant vs. variable amplitude, output options).

## Repository layout

| Path | Contents |
|---|---|
| `src/ddsim/` | The installable package — the Level I pipeline. See `docs/architecture.md` for the module map. |
| `src/ddsim/tools/` | Ported, working utilities outside the main pipeline (Weibull/reliability plotting, ANSYS output conversion). |
| `tests/` | The pytest suite, incl. `tests/fixtures/sif_verification/` (small non-proprietary verification cubes with captured original-program output). |
| `examples/` | A small runnable example case. |
| `docs/` | `architecture.md` (module map), `PORTING_NOTES.md` (the full porting history), `papers/` (reference PDFs). |
| `legacy/` | Historical, unmaintained material kept for reference: the original compiled Windows/Python-2.4 binaries (`win32-py24/`), the Windows/MPI cluster scripts being replaced by multiprocessing (`windows_cluster_scripts/`), the incomplete Level II/III prototypes (`level2_level3_prototype/`), pre-pytest ad hoc test scripts (`adhoc_test_scripts/`), and ported-but-unused modules (`unused_modules/`). |
| `reference/cpp/` | **Not tracked by git.** The original C++ Python-extension sources and their Cornell "FTools" support library, kept locally for cross-referencing while porting (see `docs/architecture.md`); they were never independently buildable (missing headers, Python 2 C-API) and have all been reimplemented in `src/ddsim`. |

## Testing

```bash
pytest
```

Two tests need real, proprietary validation data not included in this
repository (the NGC test-coupon geometry from the dissertation's validation
study) and are skipped unless `DDSIM_SIPS3002_DIR` points at a local copy —
see `docs/PORTING_NOTES.md`.

## License

MIT — see `LICENSE`.

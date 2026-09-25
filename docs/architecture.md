# Architecture

This is a map of the `ddsim` package for anyone extending it — not a tutorial on
fatigue mechanics (see the papers in `docs/papers/` and Emery et al., *Engineering
Fracture Mechanics* 76 (2009) 1500-1530, for that).

## Pipeline

```
Parameters (.par)  ─┐
mesh_io (RDB) or     ├─► MeshTools (point queries: GetPtStress,      ─┐
exodus_io (Exodus)  ─┘    IsPointOutsideMesh, surface normals, ...)   │
                                                                       ▼
                          DamMo.DamModel  ──►  DamClass (Fellipse /
                          (one instance         Hellipse / Qellipse:
                          per run, tracks        embedded / surface /
                          per-node history)      corner crack geometry,
                                                  __CalculateKi, GeometryCheck)
                                                       │
                                    GeomUtils (ellipse vs. surface-mesh    │
                                    crossings, used by GeometryCheck) ◄────┘
                                                       │
                                    Integration (RK4/RK5) + MydadN
                                    (NASGRO growth rate, Willenborg
                                    retardation)  ──►  per-node life N
                                                       │
                                    DamMo.DamModel.ToExodusFile (or
                                    ToMAPFile) ──► life as a nodal
                                    variable, e.g. for ParaView
```

`DDSim.py` (`ddsim.DDSim:main`, the `ddsim` console command) is the CLI driver:
parses arguments, builds a `MeshTools` model and a `Parameters` object, creates one
`DamMo.DamModel`, and loops over the requested node list calling either
`FwdDeterministic` (fixed initial flaw) or `MonteSimulation` (random/particle-filter
initial flaws) per node.

## Module map

| Module | Role | Tested |
|---|---|---|
| `DDSim.py` | CLI entry point / driver loop | end-to-end (`test_end_to_end.py`, `test_sif_verification_golden.py`, `test_exodus_end_to_end.py`) |
| `Parameters.py` | Reads `.par` files | yes |
| `mesh_io.py` | Reads the RDB mesh/stress ASCII format (`.con/.nod/.sig/.smp/.edg`) into `MeshData` | yes |
| `exodus_io.py` | Reads/writes Exodus II (linear elements only: `TET_4`/`WEDGE_6`/`BRICK_8`) into/from the same `MeshData`; `DamMo.ToExodusFile` writes life predictions as a nodal variable | yes, incl. against a real sample file |
| `elements.py` | Shape functions for the 6 solid + 4 surface element types (thin wrappers over `_kernels.py`) | yes |
| `_kernels.py` | numba-compiled shape functions, Newton point-location solver, bounding-box grid | yes (both directly and via `elements.py`/`MeshTools.py`) |
| `MeshTools.py` | The mesh query engine: `GetPtStress`, `IsPointOutsideMesh`, surface detection/normals, `GetNodeInfo`, etc. Replaces the compiled 2007 `MeshTools.pyd` | yes |
| `GeomUtils.py` | Ellipse-vs-surface-mesh intersections (crack front breaking the surface), `EllipseArcLength`. Replaces `GeomUtils.pyd` | yes |
| `_geom_kernels.py` | numba-compiled point-in-triangle / crossing kernel behind `GeomUtils.py` (a plain-numpy reference implementation is kept alongside it for cross-checking) | yes |
| `Vec3D.py`, `ColTensor.py`, `JohnsVectorTools.py` | Vector / symmetric-tensor / list arithmetic. Replace `Vec3D.pyd`, `ColTensor.pyd`, `JohnsVectorTools.pyd` | yes |
| `DamClass.py` | The three crack-geometry classes (`Fellipse` embedded, `Hellipse` surface, `Qellipse` corner/edge) and their stress-intensity-factor and growth logic | yes, incl. golden-output regression against the pre-port 2006 program |
| `DamMo.py` | `DamModel`: owns the mesh + per-node damage history, drives `GrowDam`/`SimDamGrowth`; `ToMAPFile`/`ToExodusFile` write results | yes |
| `DamHistory.py`, `DamErrors.py` | Small history container / exception classes | yes |
| `MydadN.py` | NASGRO crack-growth rate + Willenborg retardation (`dadN`/`Willenborg` classes). Replaces the compiled `dadN.pyd` | yes |
| `Newton.py` | Generic Newton solve used by `Willenborg`'s effective-R iteration | yes |
| `Integration.py` | Euler/RK4/RK5(Cash-Karp) steppers | yes (plus a documented anomaly in the original's error-estimate coefficients, see `docs/PORTING_NOTES.md`) |
| `Statistic.py` | Sample statistics used for Monte Carlo initial flaw sizes | yes |
| `VarAmplitude.py` | Variable-amplitude load spectrum reader | yes |
| `tools/PDF.py` | Weibull / reliability plotting (reproduces Figs. 25/27/29 of the 2009 paper) | no |
| `tools/twins.py` | ANSYS→`.sig` conversion; merges per-process output files from the old parallel workflow | no |
| `tools/rdb_to_exodus.py` | Converts an RDB mesh (+ its `.sig` stress) to Exodus for visualization -- no crack growth involved. `python -m ddsim.tools.rdb_to_exodus <rdb_base> <out.exo>` | yes |
| `tools/n_to_exodus.py` | Aggregates a real `.N` life-prediction result file (real captured 2007 Monte Carlo output, `doid rid N` per line) onto its RDB mesh as a nodal variable -- per-node arithmetic mean, matching `Statistic.StatN.SampleMean`; no crack growth run. `python -m ddsim.tools.n_to_exodus <rdb_base> <n_file> <out.exo> [var_name]` | yes |

`legacy/unused_modules/` (`Contour.py`, `Interpolator.py`, `InterpTable.py`,
`SurfaceMesh.py`) are ported (2to3'd, they compile) but are not imported by anything
above and have no test coverage. They are kept for reference, not deleted, in case
they turn out to still be wanted.

## Why the numba split (`_kernels.py` / `_geom_kernels.py`)

The mesh point-query and ellipse/surface-crossing code is the hot path (profiled:
~90% of a run's time before optimization). Rather than hand-optimize the readable
`elements.py`/`MeshTools.py`/`GeomUtils.py` code, the numeric core of each was
extracted into a separate numba-compiled module with the same algorithm, and the
readable modules became thin wrappers around it. This keeps exactly one
implementation of each algorithm (no drift between a "fast" and a "readable"
version) while the tests validate the compiled code through the same public API
used everywhere else. `GeomUtils.py` additionally keeps a pure-numpy reference
implementation (`_crossings_reference`) that the tests check against the compiled
kernel on random inputs.

## Known, intentional differences from the original (2007) code

See `docs/PORTING_NOTES.md` for the full list with rationale (surface normal
orientation, a `NearestPoint` sign bug fix, the `EllipseArcLength` step size, a
bounded rotation retry that used to recurse without limit, ...). None of them
change the golden-output regression results in `tests/test_sif_verification_golden.py`
except the one documented, understood exception (`ab3333`).

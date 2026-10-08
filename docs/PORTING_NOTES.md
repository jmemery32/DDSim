# Python 3 porting notes

Deviations from, and anomalies found in, the 2007 code. The port is meant to be
faithful: nothing below that is marked **open** has been silently changed.

## Behavior changes made during the port

| File | Change | Why |
|---|---|---|
| `VarAmplitude.py` | Load pairs are sorted **numerically** before becoming `(min, max)` | The original sorted the *strings*, so `"10.0" < "9.0"` flipped the pair. Identical results whenever the old code was right. |
| `Newton.py` | Self-test residual sign flipped (`xt - x`) | The demo diverged (`Solve` expects `F(x)` with Jacobian `dF/dx`). `Solve` itself is unchanged; only demo code. |
| `Statistic.py` | `dict.keys()/values()` wrapped in `list(...)`; sorted lists via `.sort()` on the copy | Python 3 views. Statistics unchanged; the `repr` order is still sorted. |

## Original anomalies, left as-is (**open** — decide later)

All in `Integration.py`, in the Cash–Karp RK5 error estimate. They affect
adaptive **step-size control**, not the RK5 solution itself.

1. `125.0/584` should be `125.0/594` (Cash–Karp `c4 - c4*`), in both
   `RK5_CK` and `RK5_CKslope_vector`.
2. `RK5_CKslope_vector` **adds** the `277/14336 * k5` error term; `RK5_CK`
   correctly subtracts it (`c5 = 0`, `c5* = 277/14336`).
3. `RK5_CK` uses a signed `error` (no `abs`), so a negative estimate is
   ignored by the `max(...)` floor.
4. `PEC_AM3` calls `function(t, y)` without `args`, with `t`/`y` swapped in the
   corrector evaluation. Looks unused.

## Known quirks preserved

* `Statistic.Stata.UpdateSamples` tests `if RID:`, so an explicit RID of 0
  falls into the auto-numbering branch. Harmless when RIDs start at 0.
* `Parameters` raises `IndexError` on a blank line in a `.par` file.

## Native (Python 2 C-API) modules and their replacements

| 2007 module | Replacement | Status |
|---|---|---|
| `dadN.pyd` (`dadN.dadN(mat, surf)` is really the C++ `Willenborg` class) | `MydadN.Willenborg` (`Compute_dadN` alias). `DamClass` builds these directly. | done, tested |
| `Vec3D.pyd` | `Vec3D.py` (semantics copied from `Vec3DModule.cpp`: `v*v` is dot, `Normalize` returns a new vector) | done, tested |
| `ColTensor.pyd` | `ColTensor.py` (`PrincipalValues` sorted descending, like the Jacobi routine) | done, tested |
| `JohnsVectorTools.pyd` | `JohnsVectorTools.py` (source never recovered; inferred from call sites) | done, tested |
| `MeshTools.pyd` (`MeshTools.MeshTools(name,'RDB')`) | `MeshTools.py` (queries) + `elements.py` (6 solid + 4 surface element types) + `mesh_io.py` (`MeshData`, RDB reader; Exodus will produce the same `MeshData`) | done, tested |
| `GeomUtils.pyd` (`BuildSurfMeshCObject`, `EllipseCMeshIntersections`, `EllipseArcLength`) | `GeomUtils.py` (closed-form ellipse/plane crossings, triangulated surface mesh) | done, tested |

Other things the compiled `dadN.pyd` did differently from `MydadN.py`, now aligned:
`Willenborg` substituted the material R only when the R passed in was `> 100`
(the Python version replaced any falsy R, including a legitimate `0.0`).

## Mesh layer: deliberate differences from `MeshTools.cpp`

* **Surface normals are oriented outward** using the parent element. In the
  2007 face tables the `WEDGE_6` normals point *inward* for a right-handed
  element (the other five types point outward). `DamClass` documents its
  normals as outward.
* **`NearestPoint` fixed for simplex elements** (`TET_4`, `TET_10`, `WEDGE_6`,
  `WEDGE_15`). The C++ moved the coordinates by `-u/3` with `u < 0` (wrong
  direction) and the 10-node tet mixed up `r` and `s`. Only affects the
  reported distance / element choice for points that are just outside the mesh
  (status `-1`).
* **Shape-function derivatives** come from complex-step differentiation of the
  shape functions instead of hand-typed tables (exact for polynomials).
* **Range tree -> uniform grid** over the (10%-padded) element bounding boxes.
* Not ported (unused by DDSim Level I): `GetPtDisp` (displacements were a
  hard-wired zero), `GetAdjacentSurfEdgeLengths`, `PlotSurfaceMesh`, pyramid
  elements (the C++ asserted on them in `BuildSurfaceMesh` anyway).
* `DamClass.py` (Qellipse, ~line 3871) calls `model.IsPointIn(...)`, which is not
  in the compiled module's method table -- probably dead code; to be checked.
* The repo root still contains the untracked reference directories `MeshTools/`,
  `GeomUtils/`, ... next to the new `MeshTools.py`/`GeomUtils.py` modules. Python
  prefers the `.py` file over a directory without `__init__.py`, so imports work,
  but consider moving the C++ reference sources to e.g. `legacy/cpp/`.

## Geometry layer (`GeomUtils.py`): differences from `GeomUtils.cpp`

* **Ellipse/plane crossing is closed form.** The C++ bracketed with three sample
  angles and bisected 8-14 times; it could miss thin intersections and had a
  copy-paste slip (`abs(D3) > abs(D1)` twice). Tested against dense sampling.
* **`EllipseArcLength` step size.** The 2007 default (9-point Gauss over pi/2
  steps) is off by ~0.13% for a 3:1 ellipse, more for slender ones. The step now
  scales with the aspect ratio (<1e-8 relative error up to 20:1). The integrand
  is the 2007 one. **This changes arc lengths slightly vs. 2007 results.**
* **Shared-edge crossings are reported once** (the C++ could double-count, then
  let its tangency filter delete both, or drop the real crossing).
* **The tangency/sliver filter is kept** (crossings closer than 10% of the
  perimeter along the ellipse, or wrapping >90%, are treated as a tangent touch
  and dropped, with the same "(0,1) pair => drop both, else drop first"
  rule) but its bookkeeping is now aligned; the C++ indexed `ThCross` after
  entries had been removed from `ThList`. **Heads-up:** in random tests ~45% of
  ellipses that poke through a surface have some crossings dropped by this
  heuristic -- it is a modelling choice of the original code, worth revisiting.
* `PointInFacet` decides on-edge points explicitly (the angle-sum test can cancel
  from rounding noise there); `DistanceToLineSeg` was not ported (its middle
  branch returned an uninitialised value in the C++).
* Debug `cout` output of the C++ (it printed every crossing) is gone.
* `IsPointIn` at `DamClass.py:3871` -- still to be checked.

## First end-to-end status

`examples/example1/` (12 quadratic tets, uniform syy = 3.4, `example1.par`) runs
through `DDSim.py` under Python 3.11:

    cd examples/example1
    python ../../DDSim.py -base example1 -conpath ./ -parpath ./ -v -doid_list 10 -scale 60

`tests/test_sif_validation.py` checks the stress intensity factors against
closed-form solutions: embedded elliptical cracks match Irwin to 4 decimals and
the surface half-penny matches Newman-Raju to 5 decimals (the 2007 surface
solution *is* Newman-Raju).

Other port fixes found by running it: `time.clock()` (removed in Python 3.8),
exceptions are not indexable (`message[0]` -> `message.args[0]`, 16 sites).
`Fellipse.GrowDam` in `DamClass.py` is dead code (only `DamMo.GrowDam` is called).
`Fellipse`/`Hellipse` are chosen by `GeometryCheck` (surface nodes become
`Hellipse`), not by the caller.

## Performance work (numba)

Everything below leaves the golden life predictions unchanged (bit-identical
for nodes 10 and 0, checked before and after each step; `tests/test_end_to_end.py`).

| step | example1 node 10 (scale 100) | node 0 (scale 60) |
|---|---|---|
| pure Python/numpy (before) | 7.5 s | ~139 s |
| + compiled shape functions / Newton / point query (`_kernels.py`) | 1.7 s | 5.1 s |
| + one batched stress call per SIF (`__FitPolynomial`) | 1.4 s | 3.6 s |
| + compiled ellipse/triangle crossings (`_geom_kernels.py`) | 1.2 s | 3.7 s |
| + fused normal-stress sampling (`MeshTools.NormalStressSamples`) | 1.2 s | 2.9 s |

(Wall time including ~0.5 s of interpreter/numba-cache start-up.)  All 35 nodes of
example1 at scale 100 run serially in 23.5 s.

Design: the readable numpy versions are kept next to the kernels (e.g.
`GeomUtils._crossings_reference`, `elements.py` wrappers) and tests compare the two.
`numba` caches compiled code in `__pycache__/`; the first run after a change to a
kernel takes a few extra seconds to recompile.

## Robustness fix (behavior change)

`Fellipse.GeometryCheck` recursed without bound when the whole ellipse lay outside
the body ("rotate to the second principal stress and recurse"). There are only two
orientations (sigma_1 / sigma_2 normal), so it now stops after two attempts and
returns status 4 (net fracture, crack outgrew the body) like the other "outgrew the
body" paths. Under uniaxial stress (sigma_2 = sigma_3, as in example1) the 2007 code
died with a RecursionError for any interior crack that outgrew the cube.

## Verification worth knowing about

* Cube symmetry: the eight corner nodes of example1 give the same life to ~1e-10
  relative (a direction-dependent bug would break this) -- `test_cube_symmetry_*`.
* `np.float64 * Vec3D` used to return an ndarray (numpy saw a sequence); fixed with
  `__array_ufunc__ = None` on `Vec3D` / `ColTensor`. It only bit on corner nodes.

## Golden-output verification against the original 2006 program

`tests/fixtures/sif_verification/` (copied from `SIF_Verification/{Fellipse,Hellipse}`,
small non-proprietary verification cubes, ~1000 elements) holds the literal
captured stdout of the pre-port, Windows/Python-2.4 DDSim.py run with `-verify`,
for `Fellipse` (embedded) and `Hellipse` (surface) cracks against Abaqus stress
fields of increasing polynomial order (constant, linear, bilinear, biquadratic).

`tests/test_sif_verification_golden.py` reruns these through the CURRENT
`DDSim.py` and compares every crack length and K value in the growth history,
not just the final life. **10 of 11 cases match the original to 4+ significant
figures**, including the bilinear and biquadratic (non-uniform) stress fields --
strong evidence the ported NASGRO/Willenborg/RK5/mesh-interpolation/
principal-stress chain is faithful.

### The one exception: Hellipse, depth > surface half-length ("ab3333")

The original's output for this case (`a=0.3333` depth, `b=1.0` surface
half-length) gives EXACTLY the same life as the "ab3" case (`a=1.0`, `b=0.3333`)
-- 2035.22370851 to 9 significant figures -- even though these are genuinely
different crack shapes (one wide-and-shallow along the free surface, the other
narrow-and-deep into the material). Investigated and ruled out:

* The Newman-Raju formula's two branches (`a<=c` / `a>c`) agree at the boundary
  `a==c` (continuity holds) -- `test_hellipse_newman_raju_continuous_at_depth_equals_surface_length`
  in `test_sif_validation.py`.
* The geometric `Fellipse` -> `Hellipse` conversion (`__FtoH`) does not produce
  swapped dimensions for these two inputs -- traced with instrumentation; the
  crack plane here is exactly axis-aligned (principal stress along one edge,
  surface normal along another) and there is no rotation that maps one crack
  onto the other.

Conclusion: the exact agreement in the original is most likely a genuine bug in
the 2006 code (very likely something that silently ignored which of a/b was
larger for this crack type), not a real symmetry. This port does NOT reproduce
it; `test_ab3333_known_divergence` locks in the port's own (internally
consistent) answer instead, so the difference is visible rather than silently
drifting. If this is ever confirmed against the original source/binary, update
this note.

## SIPS3002 model sanity check (real validation geometry, external data)

`tests/test_sips3002_sanity.py` loads the real SIPS3002 open-hole mesh (kept
OUTSIDE the repo -- proprietary NGC coupon geometry; set `DDSIM_SIPS3002_DIR`
to run it, otherwise skipped) through the ported `MeshTools`/`mesh_io` and
checks it against Sec. 6.1 of Emery et al. (2009):

* **140,024** `BRICK_8` + **184** `WEDGE_6` = 140,208 elements, **172,601** nodes,
  **63,974** surface nodes -- all EXACT matches to the paper.
* The 25 highest-`sigma_yy` nodes form one tight spatial cluster (<2 length
  units across, out of a model spanning ~28x6.6x1.4), consistent with the
  paper's single reported hot-spot ("aft side of hole 16 near the intersection
  of the counter-bore with the main-bore") rather than scattered noise.

Loading the real 140k-element mesh takes ~11s (one-time cost; the numba grid/
Newton kernels were only previously exercised on the 12-element `example1` and
the ~1000-element SIF_Verification cubes). This is also the node-count/surface
detection validation for the Level II/III style large models going forward.

## Real 2007 results visualized directly (`n_to_exodus.py`)

`SIPS_data/DDSimLI/SIPS3002_open/ConstantAmplitude/10000_Particles/` has real
2007 parallel-run output (`.N`/`.ai`/`.ori`) from the actual validation study.
`3002_open_CA_10000Parts_070530.N` (the latest of several dated runs in that
directory) has exactly **63,974** unique node ids -- an exact match to the
paper's surface-node count (see the sanity check above), confirming it's the
full real result set, not a partial run. `n_to_exodus.py` reads it (`doid rid
N` per line, any number of particle rows per node), reduces each node to the
plain arithmetic mean of its `N` values (matching `Statistic.StatN
.SampleMean`, the same reduction `DamModel.LifeValues`/`ToExodusFile` use),
and writes it as a `life` nodal variable via `exodus_io.write_exodus` --
nodes never seeded (not surface candidates) come back NaN. Run on the real
model: 63,974/172,601 nodes have a life value, mean 98,403, min **25,293**
cycles -- the same order of magnitude as node 113831's individual worst
particle results below, consistent (this is a mean across all particles at
each node, not the single worst particle at the single worst node, so an
exact match isn't expected). Verified end-to-end with `pvpython` (see the
Exodus section above) alongside `stress`, which `write_exodus` includes
automatically from the RDB mesh's own `.sig`.

## Planned: SIPS3002 per-particle validation (not yet run)

Node **113831** (one of the hole-16 hot-spot cluster) has 8,812 individual
particle results (life 10,904-97,112 cycles) with real `(rid, ai)` pairs
recoverable from `sips3002.map`/`sips3002.rnd` -- rich per-sample ground truth,
not just a mean. `DDSim.MonteSimulation` (already ported) can replay these
directly without the `-DB`/SQL/MPI machinery, since it only needs the ais/
ais_map file contents, not a live database. Deferred until after the
multiprocessing rewrite (§ next section) so the harness only needs to be built
once, against the final per-node execution path.

## Exodus II read/write (`exodus_io.py`)

Added ahead of multiprocessing, at the user's request, so results can be
visualized as contours in ParaView. Implemented directly against `netCDF4`
(already a dependency) rather than adding a new one; scope is linear
elements only (`TET_4`/`WEDGE_6`/`BRICK_8`), since every real DDSim mesh
(the SIPS3002 validation model, `example1`, the SIF_Verification cubes) uses
only these -- quadratic elements (`TET10`/`WEDGE15`/`HEX20`) raise a clear
`NotImplementedError` rather than risk a silently-wrong node order with no
way to verify it.

### Node-order derivation and verification

Exodus's node order does not match `elements.py`'s for `BRICK_8`: Exodus
goes CCW around one face then the opposite face in the same order; ours is
a binary-counting order (node k = `4r + 2s + t`). Rather than trust web
documentation (fetching the SEACAS/Cubit reference pages lost the node-order
diagrams through HTML->markdown summarization), the mapping was derived
empirically: a real mesh-only Exodus file was found on disk
(`Code/mpcMaker/src/TestFiles/1x2x10_beam_Left.g`, an unrelated project of
the user's, copied into `tests/fixtures/exodus/beam_hex8.g`), its
`coord`/`connect1` dumped directly, and the resulting node-to-corner mapping
hand-verified against `elements.Brick8`'s own (already-tested) natural
coordinates. `TET_4`/`WEDGE_6` use the identity mapping (the standard
Exodus/VTK/CGNS convention already matches our own structure: tet node 3 is
the apex opposite the 0-1-2 base, wedge nodes 0-1-2/3-4-5 are the bottom/top
triangles) -- no real sample file was available for these, a materially
weaker form of verification than BRICK_8's.

**A real, calibration bug found and fixed during this**: the corner-order
sanity check (`_check_nonzero_volume`, run on every element on every read)
was originally written to require a *specific sign* ("positive volume"),
calibrated so `elements.py`'s own reference node order counted as positive.
That calibration was wrong for `BRICK_8`: the real file's verified-correct
mapping legitimately has ONE natural-coordinate axis reflected relative to
`elements.py`'s reference ordering (both are valid, non-degenerate elements
-- see below), so it produced the *opposite* sign and tripped the check as
a false positive. Root cause understood in the process: DDSim's actual use
of shape functions (Newton iteration for natural coordinates, then
interpolation) only needs the element's Jacobian to be *non-singular* --
it never depends on a particular handedness/sign convention, unlike e.g. a
stiffness-matrix assembly that integrates a signed Jacobian. A pure
single-axis reflection of the natural-coordinate labeling reconstructs the
exact same physical element and gives the exact same interpolated values;
it just flips the sign of a volume formula that assumes one specific
labeling. Fixed by relaxing the check to "non-zero" (catches a genuinely
garbled/self-intersecting node order, which is the real risk) rather than
"a specific sign" -- `test_check_nonzero_volume_accepts_our_own_reference_orderings`
locks this in.

### What's tested

* `tests/fixtures/exodus/beam_hex8.g` (the real file): connectivity/coordinate
  adjacency structure (not a fixed axis assumption, which a different real
  file need not share -- see the test's docstring) and non-zero volume.
* Write -> read round trips for all three supported types (geometry, element
  connectivity, and stress, which `write_exodus` writes automatically under
  the same canonical names `read_exodus` looks for -- this is what let the
  read side get tested without a real stress-bearing sample file).
* A full DDSim crack-growth run (`test_exodus_end_to_end.py`, a single-BRICK_8
  cube with `example1`'s exact physical setup, since `example1` itself uses
  quadratic tets): real computed life -> `ToExodusFile` -> read back -> exact
  match, plus confirming un-run nodes come back as NaN (masked in ParaView),
  distinct from `DamMo.HighLife` (a node that ran and never failed).

### Three real ParaView/VTK-reader bugs found converting the actual SIPS3002 model

Round-trip tests above validate `write_exodus` against `read_exodus` -- both
ours, so a shared misunderstanding of the format wouldn't be caught. Writing
the real 172,601-node SIPS3002 model and opening it in actual ParaView
(6.2.0-RC2) surfaced three bugs invisible to that round trip, found by
reproducing "opens fine, crashes/comes back empty on Apply" directly with
`pvpython` (VTK's own `vtkIOSSReader`/`vtkExodusIIReader`, no GUI needed) so
each fix could be verified against the real reader instead of guessed at:

1. **`file_size`/`coord` mismatch.** `file_size=1` ("large model") tells
   readers to expect coordinates split into `coordx`/`coordy`/`coordz`; we
   wrote `file_size=1` but a single combined `coord` variable. VTK's reader
   errored outright on Apply: `failed to locate x nodal coordinates`
   (`RequestInformation`, which doesn't read coordinates, succeeds regardless
   -- hence "opens fine").
2. **Non-positive external IDs.** `node_num_map`/`elem_num_map` were written
   straight from DDSim's native RDB IDs, which are 0-based. VTK's IOSS reader
   rejects an ID of 0 (`node/element mapping routines detected non-positive
   global id 0`) and silently drops the *entire* element block rather than
   erroring loudly -- from the GUI this looks like "opens, then Apply gives
   an empty/broken model". Fixed by shifting each of the node and element ID
   spaces independently by the minimal constant needed to make the smallest
   ID 1, only when the source data isn't already positive (so already-1-based
   data round-trips with no offset). Connectivity is unaffected -- it already
   indexed by row position, not external ID.
3. **Combined `coord` silently zeros every nodal variable.** After fixing
   (1) by writing `file_size=0` with combined `coord`, geometry read back
   correctly but every nodal variable (`life`, `stress_xx`, ...) came back as
   all-zero -- a real quirk of this reader version specific to the
   `file_size=0`/combined-`coord` combination (confirmed by writing an
   otherwise-identical minimal file with VTK's own `vtkIOSSWriter`, which
   uses split coords + `file_size=1`, and finding *that* file's nodal
   variables read correctly). Fixed by writing split `coordx`/`coordy`/
   `coordz` + `file_size=1` throughout -- the one combination verified to get
   both geometry and nodal variables right. (`read_exodus` already handled
   both conventions; fixing it turned up a fourth, smaller bug: its split-coord
   read path returned a `numpy.ma.MaskedArray` -- a netCDF4 read artifact,
   nothing is ever actually masked in our files -- which `scipy.ConvexHull`
   flatly rejects regardless of whether anything's masked. `_read_coords` now
   returns a plain `ndarray` from both branches.)

All three are now locked in as regression tests in `test_exodus_io.py`
(`test_written_ids_are_never_non_positive`,
`test_write_leaves_already_positive_ids_unshifted`,
`test_write_uses_split_coords_matching_file_size_flag`) -- structural checks
on the written file, since the repo has no ParaView/VTK dependency to run the
real reader in CI. The regenerated `SIPS3002_open.exo` was independently
re-verified end-to-end with `pvpython` after each fix: correct point/cell
counts and the same peak `stress_xx` = 79.802928 as the earlier numeric
validation.

### Also fixed while touching `DamMo.ToMAPFile`

It referenced a bare `HighLife` instead of `self.HighLife` on its two
sentinel-value lines -- a latent `NameError` on the "no valid life computed"
path, never exercised by the test suite before now. Both `ToMAPFile` and the
new `ToExodusFile` now share one `DamModel.LifeValues(set)` helper.

### CLI

`-exodus <path>` (read mesh + stress from an Exodus file instead of RDB;
`.par` still comes from the usual `-base`/`-conpath`/`-parpath`) and
`-exodus_out <path>` (write predicted life as an Exodus nodal variable, at
the same point `-sc`/`ToMAPFile` runs). Works with `-j` too -- `-exodus`
is threaded straight through to each worker's own model construction (see
the Multiprocessing section below).

## Multiprocessing (`-j`, `parallel.py`) -- replaces the Windows/MPI cluster workflow

The original (2007) parallel execution model needed a Windows cluster: a
shared network drive, `mpirun`, an `MSTI_RANK` env var, and
`legacy/windows_cluster_scripts/*.bat` to partition a doid list into
`node_partition.<rank>` files and merge per-rank output text files back
together afterward (`ddsim.tools.twins`, still kept -- see below). `-j <N>`
replaces all of that with real multiprocessing on a single machine, in one
process's memory: partition the doid list, run each partition's crack-growth
simulation in its own worker process, merge the results, write the exact
same output files. `-p`/`-DB`/`-pn`/`-pick_list` and the code paths behind
them (`ParallelInitializeStuff`, `DataBaseInitializeStuff`,
`LocalNodeInitializeStuff`, `DoParallel`) are removed entirely --
`legacy/windows_cluster_scripts/` and `ddsim.tools.twins` are untouched
(twins is still useful for merging old already-captured cluster output, like
the real SIPS3002 `.N` files used elsewhere in this repo).

### Why a worker's results can't just cross a process boundary as-is

* `DamModel`'s per-doid container class was a **name-mangled nested class**
  (`class __DamOroContainer` inside `class DamModel`). A double-underscore
  identifier gets its *storage key* mangled by the compiler at any nesting
  depth (`DamModel.__dict__['_DamModel__DamOroContainer']`), but the class
  object's own `__qualname__` stays the unmangled `'DamModel.__DamOroContainer'`
  -- `pickle` resolves classes by walking `__qualname__` via `getattr` from
  the module, so it looks for (and doesn't find) `DamModel.__DamOroContainer`
  and raises `PicklingError`. Fixed by moving it to module level as
  `_DamOroContainer` (a single leading underscore isn't mangled at any
  nesting depth).
* Every `DamEl` entry (`DamClass.Fellipse`/`Hellipse`/`Qellipse`) holds a
  **direct reference to the entire shared `MeshTools` model**
  (`self.model = model`). Once the container rename above made a
  `DamOro[doid]` picklable *at all*, it turned out pickling one wasn't an
  *error* -- it succeeded, just catastrophically bloated: a full independent
  copy of the whole mesh gets dragged along per doid (measured: 104KB for one
  doid on a 12-element toy mesh, more than the whole mesh pickled alone,
  13.7KB -- and per-*doid*, not per-worker, on a real 172k-node mesh). But
  everything downstream only ever reads two small facts off `DamEl`:
  `WriteRotations` reads `DamEl[0].Rotation`; `WriteFinalAs` reads
  `DamEl[-1].GiveCurrent()` (an `(af,bf)` tuple). Fixed by
  `_DamOroContainer.StripDamElForTransport()`, called once a doid's
  processing is fully done (all `UpdateSample`/`CalcStats` calls complete --
  `UpdateSample` needs the live, shared `DamModel.ais` at call time to
  resolve a crack size to a particle RID, so extraction must happen after,
  not during): replaces `DamOro[doid].DamEl` with a single-element list
  holding `_StrippedDamEl`, a tiny picklable stand-in carrying just those two
  facts.
* **Accepted gaps**, both screen-diagnostics only (not file output), both
  pre-existing/documented behavior otherwise: `PrintDamInfo(dam='all')`'s
  full `DamEl` history dump comes back empty for any `-j`-processed doid;
  `FrontPoints` (already documented as "for single doid runs only") still
  runs but prints one point instead of full crack-front history for a
  `-j`-processed doid.
* `-DebugGeomUtils`/`-SVIEW` fire once, at `DamModel.__init__`, off the mesh
  itself -- not per-doid -- so the parent's own unconditionally-built
  `DamModel` in `main()` already handles them correctly regardless of `-j`;
  workers always build their own `DamModel` with these off.
* The old `errfile_name`-based error-file write (`DamModel.BuildDkva`)
  belongs exclusively to the dead `Kva_History`/`-Simp` path (disabled via
  `sys.exit()` before it would ever run). `SimDamGrowth` takes an
  `errfile_name` parameter but never uses it -- confirmed no multi-writer
  race exists in the code actually being parallelized.

### `DamMo.SimDamGrowth`'s `.Life` bug (found while designing this, fixed as a prerequisite)

`.Life` used to only be assigned inside three separate conditionals (`N >
parameters.N_max`, `N > self.HighLife`, `N < self.LowLife and will_grow !=
0`) -- so a doid whose computed life was neither capped nor a new model-wide
running max/min for the whole `DamModel` instance never got `.Life` set at
all (stayed at the container's `-1` default, later silently reported as
`self.HighLife` via `LifeValues`'s fallback, or literally `-1` via `NFile`'s
non-monte branch). Confirmed byte-for-byte identical in the very first
"Initial commit from grad school source" (`bd2d7d0`) -- an original-2007
bug, not port-introduced. Confirmed to **only affect deterministic
(non-Monte-Carlo) runs**: Monte Carlo life values come from a separate
object (`DamOro[doid].N`, via `SampleMean`) and never touch `.Life`, so none
of the already-produced/validated real SIPS3002 Monte Carlo results are
affected. Fixed by making `.Life=N` (capped the same way as before)
unconditional whenever `N>0.0`; `HighLife`/`LowLife` remain separate running
trackers. Also a correctness *prerequisite* for multiprocessing: pre-fix,
`.Life` depended on what order other doids were processed in, which doesn't
even make sense once doids are split across worker processes with their own
independent, process-local `HighLife`/`LowLife`. `DamModel.RefreshLifeBounds()`
(new) recomputes `HighLife`/`LowLife` from the merged `DamOro` after a `-j`
run, since each worker's own trackers were process-local -- cheap defensive
robustness; after this fix it's a no-op in practice, since every processed
doid's `.Life` is already correct.

### `Var_Amplitude` bugs (found and fixed alongside; NOT wired into `-j`)

Two bugs, both original-2007, both unexercised by any test before now
(`-VarAmp`/variable-amplitude loading is unvalidated/unused so far -- all
real validation is constant-amplitude):

1. References `verify` at four call sites but never receives it as a
   parameter, and no module-level `verify` global exists -- an immediate
   `NameError` on the very first `AddFDam` call, i.e. any real `-VarAmp`
   invocation crashed. Fixed by adding `verify` to its signature.
2. Its deterministic branch computed a life (`N_tot`) but never assigned
   `cracks.DamOro[doid].Life`, unlike its two Monte Carlo branches just
   above it. The naive fix (`Life = N_tot` unconditionally) is *wrong*,
   though: `VarAmp` returns `N_tot=-1` as its own legitimate sentinel for
   "reached `N_max` without growing" (setting `WillGrow` to `0`, or `-1` for
   a compressive/never-grew field) -- exactly the case the Monte Carlo
   branches already special-case (`Life = N_max` -- the local
   `1.01*parameters.N_max` -- whenever `WillGrow in (-1, 0)`, else
   `Life = N_tot`). Fixed by mirroring that same pattern in the
   deterministic branch, so `-1` from `VarAmp` never collides with
   `_DamOroContainer`'s own unrelated `-1` "never computed" sentinel.

`-VarAmp` was serial-only when this multiprocessing support was first added
(deliberately deferred at the time) -- see "`-j` wired into `-VarAmp`" further
down for when and how that changed.

### Design

`parallel.py`: `partition_doids` is the in-memory equivalent of the old
`DoParallel` -- randomly shuffles the doid list (its only load-balancing
heuristic, ported verbatim) then splits into `num_workers` contiguous,
individually-sorted chunks. A `concurrent.futures.ProcessPoolExecutor`
initializer builds each worker's own `MeshTools`+`Parameters`+`DamModel`
**once per worker process** (not once per task), amortizing the mesh-load
cost (~11s for the real SIPS3002 model) across every doid that worker
handles. Each task reuses `DDSim.MonteSimulation`/`FwdDeterministic`
directly -- the exact same dispatch `Fwd_Integration`'s serial loop body
does -- so `-j>=2` runs the identical per-doid code as `-j 1`/no `-j`, just
split across processes; `-j` absent or `1` takes today's original code path
completely unchanged. The parent merges every worker's `StripDamElForTransport`-ed
results into its own `DamModel.DamOro`, so every existing output method
(`NFile`, `WriteRotations`, `WriteInitialAs`, `WriteFinalAs`, `ToMAPFile`,
`ToExodusFile`, `LifeValues`, `NToPickle`) works unchanged against it. `-sv`
output is written once, in a batch, after all workers finish (one call per
doid, not the `'all'` batch form of `WriteInitialAs`/`WriteFinalAs`/
`WriteRotations`, which writes an extra header line the serial per-doid
loop never did) -- final file *contents* are identical to a serial `-sv`
run; the only behavior difference is that a crash mid-run no longer leaves
partial output for doids that had already completed.

### Verified

* `tests/test_parallel.py`: `partition_doids` coverage/determinism;
  `StripDamElForTransport` pickle round trip (and that pickling a live,
  unstripped `DamOro[doid]` succeeds but balloons in size -- not an error,
  a cost problem); the `SimDamGrowth`/`Var_Amplitude` bug fixes directly;
  `-j` vs. no-`-j` on `example1` (both the golden doid set and the
  8-corner-symmetry case, more doids than workers) give bit-identical
  results; `-sv` output files (`.N`/`.ai`/`.af`/`.ori`) match byte-for-byte
  between serial and `-j` runs, for both a deterministic and a Monte Carlo
  case (same `-seed`).
### Full validation against the real 2007 captured results (`test_sips3002_full_validation.py`)

An earlier draft of this section claimed a "pre-existing numerical
sensitivity" causing two independent plain serial runs of the same real
doids to disagree by several percent, based on comparing `-v` stdout prints
across two terminal invocations. **That claim was wrong** -- rechecked
properly (comparing actual `.N` file contents, not stdout) with three
independent full-precision runs on the real model's known hot-spot cluster
(two separate serial invocations, one `-j 4` run): all three produced
**bit-for-bit identical** `.N` output. There is no run-to-run
nondeterminism in the crack-growth integration, with or without `-j`.

What's actually going on, established via a full real-data comparison
(`ddsim -base SIPS3002 -conpath .../ConstantAmplitude/ -parpath
.../10000_Particles/ -S -sv -j 4`, all 63,974 real candidate nodes, compared
against the actual captured `3002_open_CA_10000Parts_070530.N`):

* **96.0%** of all 63,974 nodes match to machine precision (relative
  difference <= 1e-9); **99.4%** within 1%; **99.998%** within 10%. Mean
  *signed* relative difference across the whole model: ~1.3e-6 -- no
  systematic bias in either direction. Every doid present in one file is
  present in the other (zero missing/extra nodes both ways).
* The nodes with a meaningful difference (~350 of 63,974) are 99.4%
  concentrated in nodes with many mapped particles -- i.e. the same
  peak-stress hot-spot cluster `test_sips3002_sanity.py` already validates
  the mesh/stress reading against (Sec. 6.1). None of this ~0.6% of nodes
  falling outside 1% agreement is scattered randomly across the model.
* **Root cause, traced concretely for the worst offender** (doid 141713,
  14% off on the mean across ~6,000-8,000 particles): individual particle
  IDs (RIDs) common to both the new and the historical output do *not*
  agree well node-for-node either (e.g. RID 49839: 70,247 vs. 40,923) --
  ruling out "same inputs, different floating-point rounding" as the
  explanation. Checked why: RID 49839's initial crack size exists in the
  *current* `sips3002.rnd` but is **completely absent** from
  `sips3002.rnd.old`. Since the historical `070530.N` result reports a real
  value for this RID, **neither `.rnd`/`.map` file preserved on disk today
  is a byte-exact match to whatever snapshot actually generated that
  historical run** -- a third, no-longer-preserved version was used back in
  2007 (the `readme` in `10000_Particles/` only documents an `.old` vs.
  "current" split at "before/after May 29, 2007", not this finer-grained
  provenance gap). This also explains a separate, smaller finding: 6 of the
  63,974 nodes report "no particles nearby" (`rid=-1`, life=`N_max`) in the
  historical output while the *current* `sips3002.map` lists real particles
  for those same doids -- directly confirmed by grepping the map file.
* **Conclusion**: the rewrite (multiprocessing and the bug fixes alike) is
  correct and fully deterministic; the residual disagreement in a small,
  concentrated subset of nodes is bounded by *input-data provenance*
  (`.rnd`/`.map` files that are the best available reconstruction of 2007's
  inputs, not a verified byte-exact copy), not by anything in the code.

This is now a permanent, two-tier regression test
(`tests/test_sips3002_full_validation.py`, gated like
`test_sips3002_sanity.py` on `DDSIM_SIPS3002_DIR`): a fast (~20s) check that
actually runs `-j 2` against the real particle-filter inputs for a small
representative sample (the hot-spot cluster + 20 real sentinel nodes) and
compares against the historical result, plus a determinism check (serial vs
`-j 4` on the same real doids, asserting bit-identical output). A third,
separately-gated test (`DDSIM_SIPS3002_FULL_N`, pointing at a *precomputed*
full-model `.N` file -- it never runs the ~1hr simulation itself) checks the
whole-model statistics above against the thresholds established here.

## RK5 near-instability overshoot (2026)

Three separate, real bugs in the RK5 constant-amplitude integration path,
found while cross-checking `SimDamGrowth` (continuous RK5/Euler) against an
independent cycle-by-cycle integrator (`Var_Amplitude`) on a synthetic
near-constant-amplitude spectrum (small Gaussian noise around a fixed
amplitude, same R-ratio, on `example1`). All three are confirmed present
byte-for-byte in the original 2007 source (`bd2d7d0`), not port-introduced.

### Bug 1: Cash-Karp embedded-error coefficients (`Integration.py`)

`RK5_CK` (scalar) and `RK5_CKslope_vector` (vector) both had:

* a coefficient typo in the embedded-error combination: `125.0/584` where
  the standard Cash-Karp tableau has `125.0/594`;
* a sign error: the `k5` term in the error estimate should be *subtracted*,
  not added (`RK5_CKslope_vector` had `star(k5, 277.0/14336)`; fixed to
  `star(k5, -277.0/14336)`, matching `RK5_CK`'s own convention right next to
  it in the same file);
* a missing `abs()` before the `max()` floor check on the (signed) error
  value in `RK5_CK` -- `max(floor, error)` with a negative `error` always
  returns the floor, silently hiding real error growth.

These are real math bugs, independent of everything below, but turned out
not to be the dominant cause of the overshoot this section is mostly about
(see Bug 2) -- they're still worth having fixed regardless.

### Bug 2: the exception fallback never actually took a small step (`DamMo.GrowDam`)

`GrowDam`'s RK5 branch wraps `Integration.RK5_CKslope_vector` in a
`try/except`; on `FitPolyError`/`dAdNError`/`ValueError` it used to fall back
to `Integration.Eulerslope_vector(self.DamOro[doid].nextdN/1000.0, ...)` --
intended as "take one small step and try again next call." Two compounding
problems:

* `Eulerslope_vector`'s own `step` parameter is **completely unused in its
  body** (`return function(current_independant, current_dependant, args)`)
  -- it only evaluates the derivative at the current point, never advances
  by the given step.
* The except block only ever reassigned `self.DamOro[doid].nextdN` (pure
  bookkeeping for *next* call), never the *outer* `dN` variable that
  `SimDamGrowth` actually uses to scale the increment
  (`inc=[abs(rate[i])*dN...]`). So the "one small step" fallback silently
  used the full, unrefined `parameters.dN` as the real increment multiplier
  every time -- invisible on `example1` only because its own `dN` happens to
  equal the hardcoded `1000.0`.

Traced concretely on `example1` doid 10, scale 100: `dN` history
`[1000.0, 1.0, 1000.0, 0.0]`, crack-size history `a: 0.01 -> 0.99 -> 2.09 ->
2555.0` -- the catastrophic jump was one 1000-cycle step taken immediately
after a 1-cycle step had just shown the crack growing extremely fast. Fixed
by reassigning the outer `dN` (not just `nextdN`) to the shrunk step, and
using `parameters.dN` as the divisor instead of a hardcoded `1000.0`.

This fix alone introduced a new problem: dividing `nextdN` by
`parameters.dN` *every time* an exception recurs compounds geometrically and
can underflow to exactly `0.0`, crashing `RK5_CKslope_vector`'s `1.0/step`
with `ZeroDivisionError` (reproduced concretely on doid 8). Fixed with a
gentler, fixed divisor (`/10.0`) and a hard floor
(`parameters.dN * _MIN_STEP_FRACTION`, `_MIN_STEP_FRACTION = 0.001`).

### Bug 3 (the main fix): RK5 could take one step that overshoots the entire stable-growth region

Even with Bugs 1-2 fixed, RK5's own *successful* adaptive step-size control
could still pick a single step, near the stiff/singular region as
`Ki -> Kic`, large enough to jump the crack from "stable" to "fully
unstable, past the mesh" in one call -- not an exception path at all, just
RK5 legitimately proposing (and succeeding at) too large a step for how fast
`da/dN` is actually changing there. This is expected, textbook behavior for
explicit, fixed-order adaptive RK methods near a vertical asymptote: a local
error estimate computed for smooth local behavior doesn't anticipate
extremely fast-changing behavior just ahead of the evaluated point.

A first fix attempt capped the *starting* step proactively based on `Kmax`'s
proximity to `Kic` (a linear fraction from full step at margin >= 0.3 down to
a floor at margin -> 0). **This failed**: live-tracing doid 10 showed the
crack doubling in size in a single cycle (`a: 1.035 -> 2.013`) while this
metric still classified it "far from instability," and right at the
threshold the computed fraction was `~0.9999` -- essentially no cap exactly
where one was needed.

The working fix bounds relative crack growth *directly*, with no Ki/Kic
proxy: before each RK5 attempt, evaluate the current rate
(`DamEl.dAdN(N, abab, [scale])`) and cap the step so that

```
da/dN <= max_growth_fraction * a / dN
```

for every crack-front point with a positive rate (`parameters.max_growth_fraction`,
new `.par` keyword, default `0.1` i.e. 10% -- not named `alpha` to avoid
colliding with the existing material property `alpha`, the Willenborg
retardation shut-off exponent, `material[14]`). When this computed `safe_dN`
is smaller than the adaptive controller's current starting guess
(`self.DamOro[doid].nextdN`), RK5 is bypassed entirely in favor of a direct
Euler step at `safe_dN`. Key insight (re-confirmed by reading
`RK5_CKslope_vector`): RK5 only ever *shrinks* a step within a single call
via its own internal recursion, never grows it beyond what it was handed to
attempt -- growth only ever affects the *suggested next* call's starting
guess -- so capping the starting guess before each call is sufficient to
bound the actual step taken.

### Validation

* Independent cross-check: a synthetic near-constant-amplitude spectrum (low
  Gaussian noise around a fixed R=0 amplitude) run through both `SimDamGrowth`
  (RK5/Euler) and `Var_Amplitude` (plain cycle counting) on the same
  geometry. Pre-fix, RK5 gave life 1001.72 vs. VarAmp's 27 and plain forward
  Euler's 37.2 -- RK5 was the ~27-37x outlier. Post-fix: doid 10 -> 27.0,
  doid 8 -> 33.4, both now consistent with VarAmp/forward-Euler order of
  magnitude.
* Full real-data comparison against the historical SIPS3002 `070530.N`
  result (all 4,781 doids with real particles mapped, constant-amplitude):
  median relative difference 3.3%, p90 14.2%, max 34.4%. Variable-amplitude
  loading (403-doid representative sample) is essentially unchanged (median
  0.0%, p90 2.2%) -- expected, since `Var_Amplitude` never calls `GrowDam`
  or touches `Integration.py` at all.
* Spatial-smoothness check: the first working fix (gentle backoff divisor
  only, no growth cap) introduced visible, non-physical "splotchiness" in
  contour plots over the real SIPS3002 mesh -- confirmed quantitatively via
  a nearest-spatial-neighbor life-difference comparison on the hot-spot
  cluster (old historical result: mean local diff 764.7; that first fix:
  2854.1, ~4x worse -- genuine noise, not a visualization artifact). The
  `max_growth_fraction` cap above brings this down to 926.5 (~21% worse than
  the historical result), and the cube-symmetry test's 8 corners now agree
  to 13+ significant figures.
* `GOLDEN`/hardcoded-life values in `tests/test_end_to_end.py`,
  `tests/test_parallel.py`, and 7 of the 10 cases in
  `tests/test_sif_verification_golden.py` (moved into
  `RK5_FIX_DIVERGENCE_CASES`, same treatment as the pre-existing
  `test_ab3333_known_divergence`) were re-recorded against this fixed
  behavior -- see those files' own docstrings/comments for the reasoning per
  case. `test_cube_symmetry_...`'s 8 corners now report `WillGrow=4`, not
  the old hardcoded `2`; `test_interior_node_with_uniaxial_stress_...`
  (doid 8) now reports `WillGrow=2` instead of `4` -- the fix catches the
  same crack as genuinely unstable while it's still inside the cube, instead
  of letting it overshoot past the mesh boundary first.
* A full before/after comparison Exodus file for the real SIPS3002 model
  (`ca_life_old/new/diff`, `va_life_old/new/diff` nodal variables) was built
  and delivered to
  `SIPS_data/DDSimLI/SIPS3002_open/claude_rk5_fix_comparison/SIPS3002_old_vs_new.exo`
  for visual confirmation before these golden values were accepted.

### Still open

The ~2-3% systematic, retardation-sensitive offset found earlier in
`Var_Amplitude` validation against the historical VA result (same direction
and magnitude across the whole particle population, unlike CA's
provenance-explained scatter) was investigated (ruled out: `damp` material
parameter, both `18900` and `-nore`/no-retardation made agreement *worse*)
but not root-caused, and is independent of every fix in this section (VA
never touches `GrowDam`/RK5). Revisit if/when it matters for VA-based
results specifically.

## `-ai`/`-bi`: command-line override of the deterministic initial flaw size (2026)

Previously the only way to change a deterministic (`monte=0`) run's initial
flaw size was to hand-edit the `.par` file's `a_b` entry. `-ai <a> [-bi <b>]`
overrides it from the command line instead (`-bi` defaults to `-ai`'s value
when omitted -- the common case of an initially circular/equal-dimension
flaw), so a single geometry/`.par` pair can be swept over initial flaw size
without editing or duplicating the `.par` file. Ignored, with a printed
warning, for `monte=1`/`2` runs, where `a_b` means something different (the
sampling distribution's shape, not a single literal flaw size).

Applied once in `DDSim.main()` by overwriting `parameters.a_b` right after
the `.par` file is read -- `FwdDeterministic` (and `Var_Amplitude`'s
deterministic branch) already read the initial flaw size from
`parameters.a_b[0]`, so every call site downstream picks it up for free.
The one wrinkle: under `-j`, each worker process re-reads its own
`Parameters` from disk independently (see the Multiprocessing section
above) and would silently miss the override, re-reading the unmodified
`.par` file's own `a_b` instead -- so the override is also threaded through
as an explicit `a_b_override` argument on `parallel.run_parallel`/
`_worker_init`, applied to each worker's freshly-constructed `Parameters`
the same way the parent applies it to its own.

Verified in `tests/test_end_to_end.py`: the override actually changes the
result relative to the unmodified `.par` file; `-bi` sets the second
dimension independently; a serial run and a `-j 2` run with the same
override produce identical results; a `monte=1` run prints the warning and
otherwise runs unaffected.

## `-crack_path`: visualize the predicted crack path in ParaView (2026)

`DamModel.WriteCrackPathVTK` (wired to a new `-crack_path <path>` flag)
writes one doid's full crack-growth history as a legacy VTK PolyData file
(plain ASCII, no external dependency to write or read it -- ParaView loads
it natively) -- one polyline per recorded growth step, tracing the crack
front's actual shape in global (real, mesh) coordinates at that step, with
`CELL_DATA` giving each step's cumulative cycle count (`N`) and crack type
(`0`=Fellipse/embedded, `1`=Hellipse/surface, `2`=Qellipse/corner). Load it
in ParaView alongside the mesh/stress Exodus file (`-exodus_out`) and color
by `N` to see the predicted crack grow, step by step, over the node's life.

The per-step history and point geometry already existed (`ComputeFrontPoints`,
used by the pre-existing `-cf`/`FrontPoints`, which just dumps raw point
coordinates to a text file with no structure a viewer can use) -- the new
part is reconstructing each step's cumulative cycle count and writing it as
a structured, directly-viewable file:

* Each crack-growth regime (`DamOro[doid].DamEl`, a list -- e.g. an
  embedded Fellipse crack that later breaks the surface and becomes a
  Hellipse) keeps its own local history: `stage.a[0][0]` (a length-`L` list,
  one entry per recorded state, starting with the initial size) and
  `stage.dN` (also length `L`: `dN[i-1]` is the increment from state `i-1`
  to state `i`, for `i` in `1..L-1`; the *last* entry, `dN[L-1]`, is always
  a trailing `0.0` written by `SimDamGrowth`'s terminal/transition branches
  and never corresponds to another state in this stage -- confirmed by
  reading every code path that appends to `.dN` alongside `UpdateState`).
  Cumulative `N` is reconstructed by walking every stage in order, carrying
  the running total across stage boundaries (a transition doesn't advance
  `N` -- it's the same physical state, just re-expressed in the new
  ellipse-type's own parameterization).
* Fellipse's 20 sample points (`phi` in `[-1.0, 0.9]`, not reaching back to
  `1.0`) don't quite close the loop on their own; since Fellipse is the one
  crack type that's a fully embedded, genuinely closed ellipse (unlike
  Hellipse/Qellipse, open arcs that end at the free surface), its polyline
  repeats the first point index at the end to close it visually.

Only meaningful for one doid at a time (same "pick the first in `DamOro`"
convention as `-cf`), and only for a serial (no `-j`) run: a `-j` worker's
returned `DamOro[doid].DamEl` has already been stripped down to just the
final state (`_StrippedDamEl`/`StripDamElForTransport`, see the
Multiprocessing section above) by the time it reaches the parent -- the
full step-by-step history this needs is simply gone by then. Fails loudly
(`RuntimeError`) rather than silently writing a useless one-step file.

Verified in `tests/test_crack_path.py`: the written file's points/lines are
well-formed and self-consistent (every line indexes real points, no
duplicate indices within a line other than the deliberate Fellipse closing
repeat); the last step's cumulative `N` matches the run's own reported life;
`N` is monotonically increasing; Fellipse steps close the loop and
Hellipse/Qellipse steps don't; a `-j` run fails clearly instead of writing
a broken file.

### Bug found later: `-crack_path`'s cumulative `N` was wrong for VarAmp (2026)

The above was only tested against the constant-amplitude (`GrowDam`/RK5)
path before being shipped. Using it for real variable-amplitude SIPS3002
data (see "`ddsim.tools.top_crack_paths`" below) surfaced a real bug: the
written `N` values were off by 5-6x (e.g. a doid with a true life of 28,452
cycles showed a max recorded `N` of only 4,782).

Root cause: the original implementation reconstructed cumulative `N` by
summing `self.dN` (`cum_N += stage.dN[i-1]` for each recorded state `i`).
That's valid for `GrowDam`'s RK5/Euler path, where every `self.dN` entry but
one trailing zero lines up 1:1 with a real state -- but `VarAmp` appends to
`self.dN` on **every spectrum cycle checked**, including compressive and
near-threshold "stall" cycles that never advance `self.a[i][0]` (a crack can
stall for up to `spec.length` cycles, anywhere in the middle of a run, before
resuming growth or giving up) -- so for VarAmp, `self.dN` mixes real-growth
increments with these "nothing happened" cycles indistinguishably, and
summing it systematically undercounts.

Fixed with a new, purpose-built field: `state_N`, a list on every
`Fellipse`/`Hellipse`/`Qellipse` instance, index-aligned 1:1 with
`self.a[i][0]` (appended only on a *real* new state, in both `GrowDam` and
`VarAmp`'s growth branches -- never on a stall/skip cycle), storing the
doid's **absolute** cumulative `N` directly rather than relying on summing
anything. `WriteCrackPathVTK` now just reads `stage.state_N[i]`.

A second, related bug surfaced while fixing the first: a brand-new crack
element created at a regime transition (e.g. `Fellipse` -> `Hellipse`, same
physical state just re-expressed in the new geometry) gets a fresh
`state_N = [0.0]` from its own `__init__` -- correct for a doid's very first
element, wrong for every later one, which should inherit the current
absolute `N` at the moment of transition, not restart at zero. Fixed at both
call sites that append a new element to `DamOro[doid].DamEl` (in
`SimDamGrowth` and in `VarAmp`): immediately overwrite the new element's
`state_N` with `[N]` (the current absolute cycle count) right after
construction, before it's used.

Verified against real SIPS3002 VA data (doid 119032, see below): after both
fixes, the crack-path file's max `N` (28,451) matches the rerun's own
reported life (28,452) to within one cycle -- the one-cycle gap is
itself correct, not a bug: the final cycle that *detects* instability
(`WillGrow=2`) ends the run without growing the crack any further, so it
never gets its own recorded state.

### A real `Var_Amplitude` behavior worth knowing before using `-crack_path` on it

Discovered while building `top_crack_paths` (below): for Monte Carlo
(`monte=1/2`) variable-amplitude loading, `DDSim.Var_Amplitude` simulates
the doid's *largest* particle first (this is what sets the doid's own
reported `Life`), then walks every other particle from largest to smallest,
calling `cracks.DamOro[doid].Flush()` before each one and re-running
`VarAmp` for it -- so by the time `Var_Amplitude` returns, whatever's live
in `DamOro[doid].DamEl` is whichever particle happened to be simulated
*last* (the smallest one that still grew, right before one failed to), not
the largest/"worst case" one, and not what `Life` was actually computed
from. `-crack_path` right after a real `-VarAmp` Monte Carlo run therefore
shows an arbitrary smaller particle's path, not the doid's own headline
result -- surprising, but not a bug in `-crack_path` itself (it faithfully
reports whatever's actually left in `DamOro[doid].DamEl`). `n` other
particles in the same doid never get their own full growth history at all
(`DamModel.InterpolateLife` estimates their life from the largest
particle's curve instead) -- so there's only ever at most one real crack
path per doid.

`ddsim.tools.top_crack_paths` sidesteps this by simulating just the largest
particle directly (one `AddFDam`+`VarAmp` call, no `Var_Amplitude`/`Flush`
loop at all) -- which is also the *right* one for "most critical crack
path," since it is exactly the simulation that determines the doid's own
reported life.

## `ddsim.tools.top_crack_paths`: the N most critical crack paths (2026)

Finds the `N` lowest-life doids from an existing (real or rerun) `.N`
result file and writes each one's own representative crack-growth history
as a `-crack_path`-style VTK file -- all from a *single* mesh load, rather
than `N` separate `ddsim` subprocess invocations (each of which would pay
the real SIPS3002 mesh's own ~11s load cost independently).

For Monte Carlo (`monte=1/2`) inputs, "representative" is the doid's own
largest-particle initial flaw -- see the `Var_Amplitude`/`Flush()` behavior
above for why that specifically (not whatever `DDSim.Var_Amplitude` itself
would leave behind) is the right simulation to reproduce. Re-implements the
`.rnd`/`.map` particle-filter reading and the largest-ai `AddFDam`+`VarAmp`
(or, for constant amplitude, `SimDamGrowth`) call directly, rather than
calling `DDSim.main()`/`Var_Amplitude` -- deliberately not touching that
driver function's own, already-validated per-particle sweep.

```
python -m ddsim.tools.top_crack_paths <rdb_base> <n_file> <par_base> \
    <conpath> <parpath> <out_dir> [--val <spectrum.val>] [--top N] [--nore]
```

Verified on real SIPS3002 VA data (the historical `070530.N` as the ranking
input): the 50 lowest-life doids, reran and written as VTK crack paths,
with each rerun's own reported life matching max(`N`) in its VTK to within
one cycle (see above).

## `-j` wired into `-VarAmp` (2026)

`Var_Amplitude` finally got the same `-j <N>` multiprocessing support
`Fwd_Integration` (constant amplitude) already had.

### Design

Exactly the same shape as the constant-amplitude support (see
"Multiprocessing" above), kept as deliberately separate, parallel functions
rather than branching the existing ones:

* `DDSim.VarAmpOneDoid` -- the per-doid body factored out of
  `Var_Amplitude`'s `for doid in doid_list:` loop verbatim (same local
  variables, same logic, zero behavior change), so both the serial driver
  and the new parallel one call the *same* function -- exactly the pattern
  `Fwd_Integration` already used for `MonteSimulation`/`FwdDeterministic`.
  `Var_Amplitude` itself is now just: build the `Spectrum`, loop calling
  `VarAmpOneDoid`, do the existing `-sv` streamed-per-doid file writes --
  unchanged in substance, just a thinner wrapper.
* `parallel._worker_init_va`/`_process_chunk_va`/`run_parallel_va` --
  siblings of `_worker_init`/`_process_chunk`/`run_parallel`, not branches
  of them: the two drivers take genuinely different arguments (a spectrum
  + `nore`, no `Int_type`), and keeping them fully separate means the
  already-validated constant-amplitude path is completely untouched by
  this change. Each worker builds its own `VarAmplitude.Spectrum` in its
  initializer (cheap -- just parsing the `.val` file) alongside its own
  mesh/`Parameters`/`DamModel`. Same `StripDamElForTransport`-before-return
  memory-bounding, same merge-into-parent-`DamOro` orchestration, same
  batched `-sv` output after all workers finish.
* `DDSim.main()`'s dispatch: `-VarAmp` used to always take the serial path
  regardless of `-j` (the `elif var_file:` branch came before the
  `elif num_workers > 1:` one, so the latter was unreachable whenever the
  former matched). Now `-VarAmp` branches internally on `num_workers > 1`
  to pick `run_parallel_va` vs. the unchanged serial `Var_Amplitude`.

### Why the Monte Carlo (`monte=1/2`) branch was the real risk here

`VarAmpOneDoid`'s Monte Carlo branch mutates `DamOro[doid]` in a loop --
simulate the largest particle, then repeatedly `Flush()` and re-simulate
smaller ones until one doesn't grow (see "A real `Var_Amplitude` behavior
worth knowing" above) -- rather than computing a result in one shot the
way `MonteSimulation` does for constant amplitude. That loop had to survive
being moved, unmodified, into a function a worker process calls once per
doid in its chunk; nothing about `-j`'s chunking/merging model assumes
anything about *how* a doid's result was produced, so this worked without
needing any further change, but it was the part of this refactor worth the
most scrutiny.

### Verified

* `tests/test_parallel.py`: deterministic and Monte Carlo `-VarAmp -j 3`
  runs give bit-identical `{doid: (life, will_grow)}` vs. serial on the
  8-corner cube case (more doids than workers); `-sv` output files
  (`.N`/`.ai`/`.af`/`.ori`) match byte-for-byte between serial and `-j`
  for a deterministic run.
* Full existing suite stays green (223 passed, 5 skipped) -- in particular
  the pre-existing `-VarAmp` regression tests, confirming the
  `VarAmpOneDoid` extraction didn't change serial behavior at all.
* Against real SIPS3002 data (`monte=2`, the real particle-cracking filter,
  the trickiest case): three real hot-spot doids (119032/119036/121559,
  ~28,000 particles total across them), serial vs. `-j 3` -- the resulting
  `.N` files are byte-for-byte identical, 28,081 lines each.

## Exodus as the default output format (2026)

Previously, a run with no output flags at all wrote nothing persistent;
`-sv` opted into the old `.N`/`.ai`/`.af`/`.ori` text files, and
`-exodus_out <path>` separately opted into Exodus. Now every run writes
Exodus output by default, to `parpath+filename+'.exo'` unless
`-exodus_out <path>` picks a different path; `-sv` is unchanged (still
opts into the old text files, additively, not instead of Exodus).

Implemented as a small change to `DDSim.main()`'s existing post-dispatch
`if exodus_out:` block, which already ran after every dispatch branch
(serial or `-j`, constant or variable amplitude) once `cracks.DamOro` was
fully populated: default `exodus_out` to the auto-path when not given, and
wrap the `ToExodusFile` call in a `try/except NotImplementedError` --
*only* swallowed for the automatic default (an explicit `-exodus_out` that
fails still raises, since that was a direct request). This matters in
practice: `example1`'s own mesh uses `TET_10` (quadratic) elements, which
`exodus_io.py` doesn't support (see its own module docstring) -- without
the graceful fallback, every single existing `example1`-based use (CLI and
test) would start failing outright the moment Exodus became unconditional.
Confirmed: `example1` now prints a one-line note and continues (no file
written, run still succeeds); a real linear-element mesh (SIPS3002)
writes `<filename>.exo` and prints where. Full suite (which runs `example1`
extensively) stays green with this change -- 228 passed, 5 skipped.

### A real pre-existing bug this surfaced: the ID-shift was silently live on real data too

Found while building the companion "pull life back out of Exodus" tool
below, via a direct (non-Exodus) sanity check that didn't match: the real
SIPS3002 mesh has a node whose ID happens to be `0` (some unrelated node,
not any of the doids being checked) -- and `write_exodus`'s "shift every
ID by the minimal amount needed to make them positive" logic (added
earlier -- Exodus requires strictly positive external IDs, and ParaView/VTK
silently drops an entire element block rather than erroring on a `0`)
shifts *every* ID in the mesh together once *any* one of them is `<= 0`.
So real doids were already being silently shifted by the Exodus round trip
before this session, not just an originally-0-based toy fixture -- a doid
look-up against a `write_exodus`-produced file's raw `node_num_map` was
off by a constant for the *whole real SIPS3002 mesh*, confirmed concretely:
a direct `.N`-file mean for doid 119032 is `32844.85`, but reading it back
via the then-current `read_nodal_variable` under the key `119032` returned
`34017.33` -- the value that belonged to doid `119031`.

Fixed properly rather than just documented: `write_exodus` now records the
shift amount it used as `ddsim_node_offset`/`ddsim_elem_offset` global
netCDF attributes, and `read_exodus`/`read_nodal_variable` (both go through
the shared `_num_map` helper) undo it automatically. A file without those
attributes (written before this fix, or by something other than this
package) reads back un-adjusted, same as always -- purely additive, no
existing test needed to change (the hand-written fixtures those tests use
all already have positive IDs, so the shift -- and this fix -- is a no-op
for them; confirmed by grepping every existing use of `read_exodus`/raw
`node_num_map` reads for one that depended on the old, now-corrected,
shifted return value -- none did). Re-verified against real SIPS3002 data:
the same doid 119032 round trip now returns `32844.85`, matching the
direct `.N`-file computation exactly.

## `ddsim.tools.exodus_to_n`: the reverse of `n_to_exodus` (2026)

`n_to_exodus` goes `.N` -> Exodus; nothing went the other way. Now
`exodus_to_n` pulls one named nodal variable (default `"life"`) back out
of an Exodus file as a `.N`-style text file. Exodus stores one aggregated
value per node with no per-particle breakdown, so the output always has
exactly one line per node, `rid=-1` -- the same convention `NFile`'s own
non-Monte-Carlo branch already uses for "a single, non-particle-specific
value." A node written as NaN (`ToExodusFile`'s "never run" sentinel) is
skipped entirely, matching a doid that was never in `doid_list`, not one
that ran and got an explicit `-1`.

```
python -m ddsim.tools.exodus_to_n <exodus_file> <out.N> [var_name]
```

Built on a new, generic `exodus_io.read_nodal_variable(path, var_name)` --
unlike `read_exodus`, it doesn't need the mesh/connectivity at all, just
`node_num_map` and the matching `vals_nod_var<i>`. Building and testing this
against real data is exactly what surfaced the ID-shift bug above.

Verified: `tests/test_exodus_to_n.py` -- round-trips through `n_to_exodus`
and back (original doid numbers intact); a direct regression test for the
ID-shift fix (checks the raw on-disk IDs really are shifted, then that the
high-level reader still recovers the original ones); raises clearly on an
unknown variable name. Also checked directly against real SIPS3002 data
(see above).

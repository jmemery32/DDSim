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

## Planned: SIPS3002 per-particle validation (not yet run)

`SIPS_data/DDSimLI/SIPS3002_open/ConstantAmplitude/10000_Particles/` has real
2007 parallel-run output (`.N`/`.ai`/`.ori`) from the actual validation study.
Node **113831** (one of the hole-16 hot-spot cluster) has 8,812 individual
particle results (life 10,904-97,112 cycles) with real `(rid, ai)` pairs
recoverable from `sips3002.map`/`sips3002.rnd` -- rich per-sample ground truth,
not just a mean. `DDSim.MonteSimulation` (already ported) can replay these
directly without the `-DB`/SQL/MPI machinery, since it only needs the ais/
ais_map file contents, not a live database. Deferred until after the
multiprocessing rewrite (§ next section) so the harness only needs to be built
once, against the final per-node execution path.

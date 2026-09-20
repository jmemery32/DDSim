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
| `MeshTools.pyd` (`MeshTools.MeshTools(name,'RDB')`) | Python `Model` class + RDB reader (format documented in `MeshTools.cpp::ReadRelationalDBFiles`) | **todo** |
| `GeomUtils.pyd` (`BuildSurfMeshCObject`, `EllipseCMeshIntersections`, `EllipseArcLength`) | scipy/numpy port of `GeomUtils.cpp` | **todo** |

Other things the compiled `dadN.pyd` did differently from `MydadN.py`, now aligned:
`Willenborg` substituted the material R only when the R passed in was `> 100`
(the Python version replaced any falsy R, including a legitimate `0.0`).

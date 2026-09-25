"""numba-compiled kernels for the finite element point queries.

Everything here works on plain float64 arrays so it can be jit-compiled; the
readable, documented interface lives in :mod:`elements` and :mod:`MeshTools`,
which are thin wrappers around these functions (so there is exactly one
implementation of every shape function).

Element codes: 0 TET_4, 1 TET_10, 2 WEDGE_6, 3 WEDGE_15, 4 BRICK_8, 5 BRICK_20.

Shape functions are evaluated with *complex* arithmetic so that a complex-step
(imaginary part / h, h = 1e-30) gives exact derivatives -- see :mod:`elements`.
"""
import numpy as np
from numba import njit

TET_4, TET_10, WEDGE_6, WEDGE_15, BRICK_8, BRICK_20 = 0, 1, 2, 3, 4, 5
NUM_NODES = np.array([4, 10, 6, 15, 8, 20], dtype=np.int64)
MAX_NODES = 20
_H = 1e-30


@njit(cache=True)
def shape_c(code, r, s, t, out):
    """Shape functions at complex natural coordinates (r, s, t) -> out[:n]."""
    if code == 0:                                       # TET_4
        out[0] = r
        out[1] = s
        out[2] = t
        out[3] = 1.0 - r - s - t
    elif code == 1:                                     # TET_10
        u = 1.0 - r - s - t
        out[0] = (2 * r - 1) * r
        out[1] = (2 * s - 1) * s
        out[2] = (2 * t - 1) * t
        out[3] = (2 * u - 1) * u
        out[4] = 4 * r * s
        out[5] = 4 * s * t
        out[6] = 4 * r * t
        out[7] = 4 * t * u
        out[8] = 4 * s * u
        out[9] = 4 * r * u
    elif code == 2:                                     # WEDGE_6
        u = 1.0 - r - s
        out[0] = r * (1 - t)
        out[1] = s * (1 - t)
        out[2] = u * (1 - t)
        out[3] = r * t
        out[4] = s * t
        out[5] = u * t
    elif code == 3:                                     # WEDGE_15
        u = 1.0 - r - s
        om = 1.0 - t
        op = 1.0 + t
        t2 = 1.0 - t * t
        out[0] = 0.5 * u * ((2.0 * u - 1.0) * om - t2)
        out[1] = 0.5 * s * ((2.0 * s - 1.0) * om - t2)
        out[2] = 0.5 * r * ((2.0 * r - 1.0) * om - t2)
        out[3] = 0.5 * u * ((2.0 * u - 1.0) * op - t2)
        out[4] = 0.5 * s * ((2.0 * s - 1.0) * op - t2)
        out[5] = 0.5 * r * ((2.0 * r - 1.0) * op - t2)
        out[6] = 2.0 * u * s * om
        out[7] = 2.0 * s * r * om
        out[8] = 2.0 * r * u * om
        out[9] = 2.0 * u * s * op
        out[10] = 2.0 * s * r * op
        out[11] = 2.0 * r * u * op
        out[12] = u * t2
        out[13] = s * t2
        out[14] = r * t2
    elif code == 4:                                     # BRICK_8
        out[0] = (1.0 - r) * (1.0 - s) * (1.0 - t)
        out[1] = (1.0 - r) * (1.0 - s) * t
        out[2] = (1.0 - r) * s * (1.0 - t)
        out[3] = (1.0 - r) * s * t
        out[4] = r * (1.0 - s) * (1.0 - t)
        out[5] = r * (1.0 - s) * t
        out[6] = r * s * (1.0 - t)
        out[7] = r * s * t
    else:                                               # BRICK_20
        r2 = r * r
        s2 = s * s
        t2 = t * t
        out[0] = 0.125 * (1 - r) * (1 - s) * (1 - t) * (-r - s - t - 2.0)
        out[1] = 0.125 * (1 - r) * (1 - s) * (1 + t) * (-r - s + t - 2.0)
        out[2] = 0.125 * (1 - r) * (1 + s) * (1 - t) * (-r + s - t - 2.0)
        out[3] = 0.125 * (1 - r) * (1 + s) * (1 + t) * (-r + s + t - 2.0)
        out[4] = 0.125 * (1 + r) * (1 - s) * (1 - t) * (r - s - t - 2.0)
        out[5] = 0.125 * (1 + r) * (1 - s) * (1 + t) * (r - s + t - 2.0)
        out[6] = 0.125 * (1 + r) * (1 + s) * (1 - t) * (r + s - t - 2.0)
        out[7] = 0.125 * (1 + r) * (1 + s) * (1 + t) * (r + s + t - 2.0)
        out[8] = 0.25 * (1 - r) * (1 - s) * (1 - t2)
        out[9] = 0.25 * (1 - r) * (1 + s) * (1 - t2)
        out[10] = 0.25 * (1 + r) * (1 - s) * (1 - t2)
        out[11] = 0.25 * (1 + r) * (1 + s) * (1 - t2)
        out[12] = 0.25 * (1 - r) * (1 - s2) * (1 - t)
        out[13] = 0.25 * (1 - r) * (1 - s2) * (1 + t)
        out[14] = 0.25 * (1 + r) * (1 - s2) * (1 - t)
        out[15] = 0.25 * (1 + r) * (1 - s2) * (1 + t)
        out[16] = 0.25 * (1 - r2) * (1 - s) * (1 - t)
        out[17] = 0.25 * (1 - r2) * (1 - s) * (1 + t)
        out[18] = 0.25 * (1 - r2) * (1 + s) * (1 - t)
        out[19] = 0.25 * (1 - r2) * (1 + s) * (1 + t)


@njit(cache=True)
def shape_dshape(code, r, s, t, N, dN, tmp):
    """Real shape functions N[k] and derivatives dN[k, j] = dN_k / dnc_j.

    ``tmp`` is a complex128 scratch array of length >= 20.
    """
    n = NUM_NODES[code]
    shape_c(code, complex(r, _H), complex(s, 0.0), complex(t, 0.0), tmp)
    for k in range(n):
        N[k] = tmp[k].real
        dN[k, 0] = tmp[k].imag / _H
    shape_c(code, complex(r, 0.0), complex(s, _H), complex(t, 0.0), tmp)
    for k in range(n):
        dN[k, 1] = tmp[k].imag / _H
    shape_c(code, complex(r, 0.0), complex(s, 0.0), complex(t, _H), tmp)
    for k in range(n):
        dN[k, 2] = tmp[k].imag / _H


@njit(cache=True)
def shape_only(code, r, s, t, N, tmp):
    shape_c(code, complex(r, 0.0), complex(s, 0.0), complex(t, 0.0), tmp)
    for k in range(NUM_NODES[code]):
        N[k] = tmp[k].real


@njit(cache=True)
def center_coords(code, out):
    if code <= 1:
        out[0] = 0.25
        out[1] = 0.25
        out[2] = 0.25
    elif code == 2:
        out[0] = 1.0 / 3.0
        out[1] = 1.0 / 3.0
        out[2] = 0.5
    elif code == 3:
        out[0] = 1.0 / 3.0
        out[1] = 1.0 / 3.0
        out[2] = 0.0
    elif code == 4:
        out[0] = 0.5
        out[1] = 0.5
        out[2] = 0.5
    else:
        out[0] = 0.0
        out[1] = 0.0
        out[2] = 0.0


# ---------------------------------------------------------------------------
# inside test / nearest point
# ---------------------------------------------------------------------------
@njit(cache=True)
def point_inside(code, nc, tol):
    r, s, t = nc[0], nc[1], nc[2]
    if code <= 1:                                       # tetrahedra
        u = 1.0 - r - s - t
        return (r >= -tol and r <= 1.0 + tol and s >= -tol and s <= 1.0 + tol and
                t >= -tol and t <= 1.0 + tol and u >= -tol and u <= 1.0 + tol)
    if code == 2 or code == 3:                          # wedges
        u = 1.0 - r - s
        ok = (r >= -tol and r <= 1.0 + tol and s >= -tol and s <= 1.0 + tol and
              u >= -tol and u <= 1.0 + tol)
        if code == 2:
            return ok and t >= -tol and t <= 1.0 + tol
        return ok and t >= -(1.0 + tol) and t <= 1.0 + tol
    if code == 4:                                       # BRICK_8, [0, 1]^3
        return (r >= -tol and r <= 1.0 + tol and s >= -tol and s <= 1.0 + tol and
                t >= -tol and t <= 1.0 + tol)
    lim = 1.0 + tol                                     # BRICK_20, [-1, 1]^3
    return (r >= -lim and r <= lim and s >= -lim and s <= lim and t >= -lim and t <= lim)


@njit(cache=True)
def _project_simplex(c, n, out):
    """Euclidean projection of c[:n] onto {x >= 0, sum(x) <= 1}."""
    total = 0.0
    for i in range(n):
        out[i] = max(c[i], 0.0)
        total += out[i]
    if total <= 1.0:
        return
    u = np.sort(out[:n])[::-1]
    css = 0.0
    theta = 0.0
    for k in range(n):
        css += u[k]
        cand = (css - 1.0) / (k + 1.0)
        if u[k] - cand > 0.0:
            theta = cand
    for i in range(n):
        out[i] = max(out[i] - theta, 0.0)


@njit(cache=True)
def nearest_point(code, nc, out):
    """Nearest natural-coordinate point inside the reference element."""
    if code <= 1:
        _project_simplex(nc, 3, out)
    elif code == 2 or code == 3:
        _project_simplex(nc, 2, out)
        out[2] = min(max(nc[2], 0.0 if code == 2 else -1.0), 1.0)
    elif code == 4:
        for i in range(3):
            out[i] = min(max(nc[i], 0.0), 1.0)
    else:
        for i in range(3):
            out[i] = min(max(nc[i], -1.0), 1.0)


# ---------------------------------------------------------------------------
# point location
# ---------------------------------------------------------------------------
@njit(cache=True)
def find_natural(code, X, q, tol, max_its, nc, N, dN, tmp):
    """Newton iteration for the natural coordinates ``nc`` of point ``q``.

    ``X`` is the (n, 3) array of the element's node coordinates.  Returns True
    on convergence (``nc`` holds the result), False otherwise.
    """
    n = NUM_NODES[code]
    center_coords(code, nc)
    for _ in range(max_its):
        shape_dshape(code, nc[0], nc[1], nc[2], N, dN, tmp)
        g0 = 0.0
        g1 = 0.0
        g2 = 0.0
        j00 = j01 = j02 = j10 = j11 = j12 = j20 = j21 = j22 = 0.0
        for k in range(n):
            x, y, z = X[k, 0], X[k, 1], X[k, 2]
            g0 += N[k] * x
            g1 += N[k] * y
            g2 += N[k] * z
            j00 += dN[k, 0] * x
            j01 += dN[k, 0] * y
            j02 += dN[k, 0] * z
            j10 += dN[k, 1] * x
            j11 += dN[k, 1] * y
            j12 += dN[k, 1] * z
            j20 += dN[k, 2] * x
            j21 += dN[k, 2] * y
            j22 += dN[k, 2] * z
        d0 = g0 - q[0]
        d1 = g1 - q[1]
        d2 = g2 - q[2]
        # jac[i][j] = d x_j / d nc_i ;  solve jac^T * update = delta  (Cramer's rule)
        det = (j00 * (j11 * j22 - j21 * j12) - j10 * (j01 * j22 - j21 * j02) +
               j20 * (j01 * j12 - j11 * j02))
        if det == 0.0 or not np.isfinite(det):
            return False
        u0 = (d0 * (j11 * j22 - j21 * j12) - j10 * (d1 * j22 - j21 * d2) +
              j20 * (d1 * j12 - j11 * d2)) / det
        u1 = (j00 * (d1 * j22 - j21 * d2) - d0 * (j01 * j22 - j21 * j02) +
              j20 * (j01 * d2 - d1 * j02)) / det
        u2 = (j00 * (j11 * d2 - d1 * j12) - j10 * (j01 * d2 - d1 * j02) +
              d0 * (j01 * j12 - j11 * j02)) / det
        n0, n1, n2 = nc[0] - u0, nc[1] - u1, nc[2] - u2
        if not (np.isfinite(n0) and np.isfinite(n1) and np.isfinite(n2)):
            return False
        moved = abs(n0 - nc[0]) > tol or abs(n1 - nc[1]) > tol or abs(n2 - nc[2]) > tol
        nc[0], nc[1], nc[2] = n0, n1, n2
        if not moved:
            return True
    return False


@njit(cache=True)
def locate(code, conn_row, xyz, q, tol, max_its, nc, X, N, dN, tmp):
    """Gather element node coordinates and run :func:`find_natural`."""
    n = NUM_NODES[code]
    for k in range(n):
        row = conn_row[k]
        X[k, 0] = xyz[row, 0]
        X[k, 1] = xyz[row, 1]
        X[k, 2] = xyz[row, 2]
    return find_natural(code, X, q, tol, max_its, nc, N, dN, tmp)


@njit(cache=True)
def interpolate(code, conn_row, field, nc, out, N, tmp):
    """out[c] = sum_k N_k(nc) * field[conn_row[k], c]   (field is (n_nodes, ncomp))."""
    shape_only(code, nc[0], nc[1], nc[2], N, tmp)
    ncomp = field.shape[1]
    for c in range(ncomp):
        out[c] = 0.0
    for k in range(NUM_NODES[code]):
        row = conn_row[k]
        w = N[k]
        for c in range(ncomp):
            out[c] += w * field[row, c]


@njit(cache=True)
def surface_distance(code, conn_row, xyz, q, nc, near_nc, X, N, tmp):
    """Distance from ``q`` to the reference-element point nearest to ``nc`` mapped to space."""
    nearest_point(code, nc, near_nc)
    shape_only(code, near_nc[0], near_nc[1], near_nc[2], N, tmp)
    px = 0.0
    py = 0.0
    pz = 0.0
    for k in range(NUM_NODES[code]):
        row = conn_row[k]
        px += N[k] * xyz[row, 0]
        py += N[k] * xyz[row, 1]
        pz += N[k] * xyz[row, 2]
    return np.sqrt((q[0] - px) ** 2 + (q[1] - py) ** 2 + (q[2] - pz) ** 2)


# ---------------------------------------------------------------------------
# uniform grid over the element bounding boxes (CSR layout)
# ---------------------------------------------------------------------------
@njit(cache=True)
def build_grid(lo, hi, origin, inv_cell, dims):
    """Register every element in each grid cell its bounding box touches.

    Returns ``(cell_start, cell_items)``: the elements of cell ``c`` are
    ``cell_items[cell_start[c]:cell_start[c + 1]]`` in ascending element order.
    """
    m = lo.shape[0]
    ncell = dims[0] * dims[1] * dims[2]
    counts = np.zeros(ncell + 1, dtype=np.int64)
    for pas in range(2):
        if pas == 1:
            start = np.zeros(ncell + 1, dtype=np.int64)
            for c in range(ncell):
                start[c + 1] = start[c] + counts[c]
            items = np.empty(start[ncell], dtype=np.int64)
            fill = start[:ncell].copy()
        for i in range(m):
            ilo0 = int(np.floor((lo[i, 0] - origin[0]) * inv_cell))
            ilo1 = int(np.floor((lo[i, 1] - origin[1]) * inv_cell))
            ilo2 = int(np.floor((lo[i, 2] - origin[2]) * inv_cell))
            ihi0 = int(np.floor((hi[i, 0] - origin[0]) * inv_cell))
            ihi1 = int(np.floor((hi[i, 1] - origin[1]) * inv_cell))
            ihi2 = int(np.floor((hi[i, 2] - origin[2]) * inv_cell))
            for a in range(ilo0, ihi0 + 1):
                for b in range(ilo1, ihi1 + 1):
                    for c in range(ilo2, ihi2 + 1):
                        cell = (a * dims[1] + b) * dims[2] + c
                        if pas == 0:
                            counts[cell] += 1
                        else:
                            items[fill[cell]] = i
                            fill[cell] += 1
    return start, items


# ---------------------------------------------------------------------------
# the point query
# ---------------------------------------------------------------------------
@njit(cache=True)
def pt_query(q, cached, tol, max_its, origin, inv_cell, dims, cell_start, cell_items,
             lo, hi, code, conn, xyz, nat, nc, X, N, dN, tmp, near):
    """Locate ``q`` in the mesh.  Returns ``(status, element index, distance)``.

    status  -2  inside element ``i`` (``nat`` = its natural coordinates)
            -1  just outside; ``i`` is the nearest element, ``nat`` the nearest
                point of it, ``distance`` the distance to that point
            -3  no element bounding box contains ``q``
            -4  Newton iteration did not converge for any candidate; ``i`` is
                the last such element
    ``cached`` is the element the previous query landed in (or -1).
    """
    if cached >= 0:
        if (locate(code[cached], conn[cached], xyz, q, tol, max_its, nc, X, N, dN, tmp)
                and point_inside(code[cached], nc, -tol)):
            nat[0], nat[1], nat[2] = nc[0], nc[1], nc[2]
            return -2, cached, 0.0

    ix = int(np.floor((q[0] - origin[0]) * inv_cell))
    iy = int(np.floor((q[1] - origin[1]) * inv_cell))
    iz = int(np.floor((q[2] - origin[2]) * inv_cell))
    if ix < 0 or iy < 0 or iz < 0 or ix >= dims[0] or iy >= dims[1] or iz >= dims[2]:
        return -3, -1, 0.0
    cell = (ix * dims[1] + iy) * dims[2] + iz

    any_candidate = False
    last_unconverged = -1
    best_i = -1
    best_d = 0.0
    for p in range(cell_start[cell], cell_start[cell + 1]):
        i = cell_items[p]
        if (q[0] < lo[i, 0] or q[0] > hi[i, 0] or q[1] < lo[i, 1] or q[1] > hi[i, 1]
                or q[2] < lo[i, 2] or q[2] > hi[i, 2]):
            continue
        any_candidate = True
        c = code[i]
        if not locate(c, conn[i], xyz, q, tol, max_its, nc, X, N, dN, tmp):
            last_unconverged = i
            continue
        if point_inside(c, nc, tol):
            nat[0], nat[1], nat[2] = nc[0], nc[1], nc[2]
            return -2, i, 0.0
        dist = surface_distance(c, conn[i], xyz, q, nc, near, X, N, tmp)
        if best_i < 0 or dist < best_d:
            best_i = i
            best_d = dist
            nat[0], nat[1], nat[2] = near[0], near[1], near[2]

    if not any_candidate:
        return -3, -1, 0.0
    if best_i < 0:
        return -4, last_unconverged, 0.0
    return -1, best_i, best_d


@njit(cache=True)
def stress_batch(pts, cached, tol, max_its, origin, inv_cell, dims, cell_start, cell_items,
                 lo, hi, code, conn, xyz, sig, out):
    """Interpolated stress at every row of ``pts`` (n, 3) -> ``out`` (n, ncomp).

    Stops at the first point that is not in (or just next to) the mesh.
    Returns ``(n_done, cached)``: ``n_done < n`` means point ``n_done`` failed.
    """
    nat = np.empty(3)
    nc = np.empty(3)
    near = np.empty(3)
    X = np.empty((MAX_NODES, 3))
    N = np.empty(MAX_NODES)
    dN = np.empty((MAX_NODES, 3))
    tmp = np.empty(MAX_NODES, dtype=np.complex128)
    for k in range(pts.shape[0]):
        status, i, dist = pt_query(pts[k], cached, tol, max_its, origin, inv_cell, dims,
                                   cell_start, cell_items, lo, hi, code, conn, xyz,
                                   nat, nc, X, N, dN, tmp, near)
        if status == -2:
            cached = i
        elif status == -1:
            cached = -1
        else:
            return k, cached
        interpolate(code[i], conn[i], sig, nat, out[k], N, tmp)
    return pts.shape[0], cached


@njit(cache=True)
def normal_stress_samples(R, center, ypnts, zpnts, cached, tol, max_its, origin, inv_cell,
                          dims, cell_start, cell_items, lo, hi, code, conn, xyz, sig, f):
    """Normal stress on a crack plane sampled on a grid, then at its center.

    Points are ``center + y * R[1] + z * R[2]`` for every ``(y, z)`` pair (``y``
    varying slowest) followed by ``center`` itself; ``R``'s rows are the crack
    frame (``R[0]`` = crack-plane normal).  ``f[k] = R0 . sigma(x_k) . R0``,
    evaluated with the same expression and operation order as the 2007 loops
    (``sig`` columns: xx yy zz xy yz zx).

    Returns ``(n_done, cached)``; ``n_done`` < number of points means a point
    fell outside the mesh.
    """
    ny = ypnts.shape[0]
    nz = zpnts.shape[0]
    npts = ny * nz + 1
    pts = np.empty((npts, 3))
    idx = 0
    for a in range(ny):
        for b in range(nz):
            for k in range(3):
                pts[idx, k] = center[k] + (R[1, k] * ypnts[a] + R[2, k] * zpnts[b])
            idx += 1
    for k in range(3):
        pts[idx, k] = center[k]
    out = np.empty((npts, 6))
    done, cached = stress_batch(pts, cached, tol, max_its, origin, inv_cell, dims, cell_start,
                                cell_items, lo, hi, code, conn, xyz, sig, out)
    if done < npts:
        return done, cached
    A = R[0, 0]
    B = R[0, 1]
    C = R[0, 2]
    for i in range(npts):
        f[i] = (A * A * out[i, 0] + B * B * out[i, 1] + C * C * out[i, 2]
                + 2.0 * A * B * out[i, 3] + 2.0 * A * C * out[i, 5] + 2.0 * B * C * out[i, 4])
    return npts, cached

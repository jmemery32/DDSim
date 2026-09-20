"""numba kernels for GeomUtils: ellipse vs. triangle crossings.

Same algorithms as the readable numpy versions in :mod:`GeomUtils`
(``PointInFacet``, ``_plane_crossings``, ``_ellipse_point``), specialised to the
triangles the surface mesh is made of.  ``tests/test_geomutils.py`` checks the
two implementations against each other.
"""
import math

import numpy as np
from numba import njit

_TWO_PI = 2.0 * math.pi


@njit(cache=True)
def point_in_triangle(p, tri):
    """PointInFacet for a 3-corner facet (tol = 0)."""
    f = np.empty((3, 3))
    for i in range(3):
        for k in range(3):
            f[i, k] = tri[i, k] - p[k]
    ma = np.empty(3)
    for i in range(3):
        ma[i] = math.sqrt(f[i, 0] ** 2 + f[i, 1] ** 2 + f[i, 2] ** 2)
        if ma[i] == 0.0:                                 # the point is a vertex
            return True

    # on an edge?  (the angle-sum below is noise there)
    for i in range(3):
        j = (i + 1) % 3
        sx = f[j, 0] - f[i, 0]
        sy = f[j, 1] - f[i, 1]
        sz = f[j, 2] - f[i, 2]
        seg2 = sx * sx + sy * sy + sz * sz
        if seg2 == 0.0:
            continue
        t = -(f[i, 0] * sx + f[i, 1] * sy + f[i, 2] * sz) / seg2
        t = min(max(t, 0.0), 1.0)
        cx = f[i, 0] + t * sx
        cy = f[i, 1] + t * sy
        cz = f[i, 2] + t * sz
        if math.sqrt(cx * cx + cy * cy + cz * cz) <= 1e-10 * math.sqrt(seg2):
            return True

    # winding-angle sum
    e1x, e1y, e1z = f[1, 0] - f[0, 0], f[1, 1] - f[0, 1], f[1, 2] - f[0, 2]
    e2x, e2y, e2z = f[2, 0] - f[1, 0], f[2, 1] - f[1, 1], f[2, 2] - f[1, 2]
    nx = e1y * e2z - e1z * e2y
    ny = e1z * e2x - e1x * e2z
    nz = e1x * e2y - e1y * e2x
    total = 0.0
    for i in range(3):
        j = (i + 1) % 3
        dot = f[i, 0] * f[j, 0] + f[i, 1] * f[j, 1] + f[i, 2] * f[j, 2]
        cosang = min(max(dot / (ma[i] * ma[j]), -1.0), 1.0)
        cx = f[i, 1] * f[j, 2] - f[i, 2] * f[j, 1]
        cy = f[i, 2] * f[j, 0] - f[i, 0] * f[j, 2]
        cz = f[i, 0] * f[j, 1] - f[i, 1] * f[j, 0]
        d = nx * cx + ny * cy + nz * cz
        sign = 1.0 if d > 0.0 else (-1.0 if d < 0.0 else 0.0)
        total += sign * math.acos(cosang)
    return abs(total) > math.pi / 2.0


@njit(cache=True)
def _polar_angle(a, b, psi):
    return math.atan2(b * math.sin(psi), a * math.cos(psi)) % _TWO_PI


@njit(cache=True)
def crossings(tri, normal, offset, ids, center, e1, e2, a, b, th_out, j_out, t_out):
    """All (theta, j, triangle) crossings of the ellipse with the triangles ``ids``.

    ``j`` is 0 for the smaller and 1 for the larger of a triangle's two plane
    crossings.  Results are written to the three ``*_out`` arrays (length >= 2 *
    len(ids)); the number of crossings is returned.
    """
    n = 0
    for idx in range(len(ids)):
        t = ids[idx]
        p = normal[t, 0] * e1[0] + normal[t, 1] * e1[1] + normal[t, 2] * e1[2]
        q = normal[t, 0] * e2[0] + normal[t, 1] * e2[1] + normal[t, 2] * e2[2]
        d0 = (normal[t, 0] * center[0] + normal[t, 1] * center[1]
              + normal[t, 2] * center[2]) + offset[t]
        amp = math.hypot(p * a, q * b)
        if amp == 0.0:
            continue                                     # plane parallel to the ellipse plane
        c = -d0 / amp
        if abs(c) > 1.0:
            if abs(c) > 1.0 + 1e-12:
                continue
            c = 1.0 if c > 0.0 else -1.0
        psi0 = math.atan2(q * b, p * a)
        alpha = math.acos(c)
        th0 = _polar_angle(a, b, psi0 - alpha)
        th1 = _polar_angle(a, b, psi0 + alpha)
        if th0 > th1:
            th0, th1 = th1, th0
        pt = np.empty(3)
        for j in range(2):
            theta = th0 if j == 0 else th1
            cs, sn = math.cos(theta), math.sin(theta)
            r = a * b / math.sqrt(b * b * cs * cs + a * a * sn * sn)
            for k in range(3):
                pt[k] = center[k] + r * (cs * e1[k] + sn * e2[k])
            if point_in_triangle(pt, tri[t]):
                th_out[n] = theta
                j_out[n] = j
                t_out[n] = t
                n += 1
    return n

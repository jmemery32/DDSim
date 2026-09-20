"""Finite element shape functions and reference-element geometry.

Python port of the element classes in ``ElemClass.cpp`` (MeshTools, 2007) for
the six solid element types DDSim handles, plus the four surface (face)
element types used to compute surface normals.

Node numbering, natural-coordinate conventions and face tables are identical
to the C++ code, so ``.con``/``.edg`` files and Exodus-style connectivity can be
mapped onto them unchanged:

============  =====  ===============================  ==========
element       nodes  natural coordinates              center
============  =====  ===============================  ==========
Tet4/Tet10    4/10   r,s,t >= 0, r+s+t <= 1           (1/4,1/4,1/4)
Wedge6        6      (r,s) triangle, t in [0, 1]      (1/3,1/3,1/2)
Wedge15       15     (r,s) triangle, t in [-1, 1]     (1/3,1/3,0)
Brick8        8      r,s,t in [0, 1]                  (1/2,1/2,1/2)
Brick20       20     r,s,t in [-1, 1]                 (0,0,0)
============  =====  ===============================  ==========

Shape-function derivatives are obtained by complex-step differentiation of the
shape functions.  The shape functions are polynomials, so this is exact to
machine precision and removes ~300 lines of hand-derived derivative tables (and
their chance of typos).  Every shape function below therefore uses only
arithmetic that works on complex numbers.
"""
import numpy as np

_H = 1e-30  # complex-step size; exact for polynomials


# ---------------------------------------------------------------------------
# solid elements
# ---------------------------------------------------------------------------
class SolidElement:
    """Reference solid element.  Subclasses implement ``_shape``."""

    name = ""
    num_nodes = 0
    center = (0.0, 0.0, 0.0)
    faces = ()          # tuple of node-index tuples (corner nodes first, then mid-side)

    def shape(self, nc):
        """Shape function values at natural coordinates ``nc`` -> ndarray(num_nodes)."""
        return self._shape(float(nc[0]), float(nc[1]), float(nc[2]))

    def dshape(self, nc):
        """Shape derivatives d N_k / d nc_j -> ndarray(num_nodes, 3)."""
        d = np.empty((self.num_nodes, 3))
        r, s, t = float(nc[0]), float(nc[1]), float(nc[2])
        d[:, 0] = self._shape(complex(r, _H), s, t).imag / _H
        d[:, 1] = self._shape(r, complex(s, _H), t).imag / _H
        d[:, 2] = self._shape(r, s, complex(t, _H)).imag / _H
        return d

    def point_inside(self, nc, tol):
        raise NotImplementedError

    def nearest_point(self, nc):
        raise NotImplementedError

    def _shape(self, r, s, t):
        raise NotImplementedError


def _project_onto_simplex(c):
    """Nearest point to ``c`` on {x >= 0, sum(x) <= 1} (Euclidean projection).

    The 2007 code moved the coordinates by ``-u/3`` with ``u = 1 - sum < 0`` (i.e.
    the wrong way) and, for the 10-node tet, mixed up r and s; this is the
    correct projection.  It only affects the reported distance for points that
    are (slightly) outside an element.
    """
    c = np.clip(np.asarray(c, dtype=float), 0.0, None)
    if c.sum() <= 1.0:
        return c
    # projection onto the probability simplex (sort-based algorithm)
    u = np.sort(c)[::-1]
    css = np.cumsum(u) - 1.0
    k = np.nonzero(u - css / (np.arange(len(u)) + 1) > 0)[0][-1]
    theta = css[k] / (k + 1.0)
    return np.clip(c - theta, 0.0, None)


class Tet4(SolidElement):
    name = "TET_4"
    num_nodes = 4
    center = (0.25, 0.25, 0.25)
    faces = ((0, 1, 2), (1, 3, 2), (2, 3, 0), (0, 3, 1))

    def _shape(self, r, s, t):
        return np.array([r, s, t, 1.0 - r - s - t])

    def point_inside(self, nc, tol):
        r, s, t = nc
        u = 1.0 - r - s - t
        return all(-tol <= v <= 1.0 + tol for v in (r, s, t, u))

    def nearest_point(self, nc):
        return _project_onto_simplex(nc)


class Tet10(SolidElement):
    name = "TET_10"
    num_nodes = 10
    center = (0.25, 0.25, 0.25)
    faces = ((0, 1, 2, 4, 5, 6), (1, 3, 2, 8, 7, 5),
             (2, 3, 0, 7, 9, 6), (0, 3, 1, 9, 8, 4))
    point_inside = Tet4.point_inside
    nearest_point = Tet4.nearest_point

    def _shape(self, r, s, t):
        u = 1.0 - r - s - t
        return np.array([(2 * r - 1) * r, (2 * s - 1) * s, (2 * t - 1) * t, (2 * u - 1) * u,
                         4 * r * s, 4 * s * t, 4 * r * t, 4 * t * u, 4 * s * u, 4 * r * u])


class Wedge6(SolidElement):
    name = "WEDGE_6"
    num_nodes = 6
    center = (1.0 / 3.0, 1.0 / 3.0, 0.5)
    faces = ((0, 3, 4, 1), (2, 5, 3, 0), (1, 4, 5, 2), (0, 1, 2), (3, 5, 4))

    def _shape(self, r, s, t):
        u = 1.0 - r - s
        return np.array([r * (1 - t), s * (1 - t), u * (1 - t), r * t, s * t, u * t])

    def point_inside(self, nc, tol):
        r, s, t = nc
        u = 1.0 - r - s
        return all(-tol <= v <= 1.0 + tol for v in (r, s, t, u))

    def nearest_point(self, nc):
        rs = _project_onto_simplex(nc[:2])
        return np.array([rs[0], rs[1], min(max(nc[2], 0.0), 1.0)])


class Wedge15(SolidElement):
    name = "WEDGE_15"
    num_nodes = 15
    center = (1.0 / 3.0, 1.0 / 3.0, 0.0)
    faces = ((0, 3, 4, 1, 12, 9, 13, 6), (2, 5, 3, 0, 14, 11, 12, 8),
             (1, 4, 5, 2, 13, 10, 14, 7), (0, 1, 2, 6, 7, 8), (3, 5, 4, 11, 10, 9))

    def _shape(self, r, s, t):
        u = 1.0 - r - s
        om, op, t2 = 1.0 - t, 1.0 + t, 1.0 - t * t
        return np.array([
            0.5 * u * ((2.0 * u - 1.0) * om - t2),
            0.5 * s * ((2.0 * s - 1.0) * om - t2),
            0.5 * r * ((2.0 * r - 1.0) * om - t2),
            0.5 * u * ((2.0 * u - 1.0) * op - t2),
            0.5 * s * ((2.0 * s - 1.0) * op - t2),
            0.5 * r * ((2.0 * r - 1.0) * op - t2),
            2.0 * u * s * om, 2.0 * s * r * om, 2.0 * r * u * om,
            2.0 * u * s * op, 2.0 * s * r * op, 2.0 * r * u * op,
            u * t2, s * t2, r * t2])

    def point_inside(self, nc, tol):
        r, s, t = nc
        u = 1.0 - r - s
        return all(-tol <= v <= 1.0 + tol for v in (r, s, u)) and -(1.0 + tol) <= t <= 1.0 + tol

    def nearest_point(self, nc):
        rs = _project_onto_simplex(nc[:2])
        return np.array([rs[0], rs[1], min(max(nc[2], -1.0), 1.0)])


class Brick8(SolidElement):
    name = "BRICK_8"
    num_nodes = 8
    center = (0.5, 0.5, 0.5)
    faces = ((0, 1, 3, 2), (4, 6, 7, 5), (0, 4, 5, 1),
             (2, 3, 7, 6), (0, 2, 6, 4), (1, 5, 7, 3))

    def _shape(self, r, s, t):
        return np.array([(1.0 - r) * (1.0 - s) * (1.0 - t), (1.0 - r) * (1.0 - s) * t,
                         (1.0 - r) * s * (1.0 - t), (1.0 - r) * s * t,
                         r * (1.0 - s) * (1.0 - t), r * (1.0 - s) * t,
                         r * s * (1.0 - t), r * s * t])

    def point_inside(self, nc, tol):
        return all(-tol <= v <= 1.0 + tol for v in nc)

    def nearest_point(self, nc):
        return np.clip(np.asarray(nc, dtype=float), 0.0, 1.0)


class Brick20(SolidElement):
    name = "BRICK_20"
    num_nodes = 20
    center = (0.0, 0.0, 0.0)
    faces = ((0, 1, 3, 2, 8, 13, 9, 12), (4, 6, 7, 5, 14, 11, 15, 10),
             (0, 4, 5, 1, 16, 10, 17, 8), (2, 3, 7, 6, 9, 19, 11, 18),
             (0, 2, 6, 4, 12, 18, 14, 16), (1, 5, 7, 3, 17, 15, 19, 13))

    def _shape(self, r, s, t):
        r2, s2, t2 = r * r, s * s, t * t
        return np.array([
            0.125 * (1 - r) * (1 - s) * (1 - t) * (-r - s - t - 2.0),
            0.125 * (1 - r) * (1 - s) * (1 + t) * (-r - s + t - 2.0),
            0.125 * (1 - r) * (1 + s) * (1 - t) * (-r + s - t - 2.0),
            0.125 * (1 - r) * (1 + s) * (1 + t) * (-r + s + t - 2.0),
            0.125 * (1 + r) * (1 - s) * (1 - t) * (r - s - t - 2.0),
            0.125 * (1 + r) * (1 - s) * (1 + t) * (r - s + t - 2.0),
            0.125 * (1 + r) * (1 + s) * (1 - t) * (r + s - t - 2.0),
            0.125 * (1 + r) * (1 + s) * (1 + t) * (r + s + t - 2.0),
            0.25 * (1 - r) * (1 - s) * (1 - t2), 0.25 * (1 - r) * (1 + s) * (1 - t2),
            0.25 * (1 + r) * (1 - s) * (1 - t2), 0.25 * (1 + r) * (1 + s) * (1 - t2),
            0.25 * (1 - r) * (1 - s2) * (1 - t), 0.25 * (1 - r) * (1 - s2) * (1 + t),
            0.25 * (1 + r) * (1 - s2) * (1 - t), 0.25 * (1 + r) * (1 - s2) * (1 + t),
            0.25 * (1 - r2) * (1 - s) * (1 - t), 0.25 * (1 - r2) * (1 - s) * (1 + t),
            0.25 * (1 - r2) * (1 + s) * (1 - t), 0.25 * (1 - r2) * (1 + s) * (1 + t)])

    def point_inside(self, nc, tol):
        return all(-(1.0 + tol) <= v <= 1.0 + tol for v in nc)

    def nearest_point(self, nc):
        return np.clip(np.asarray(nc, dtype=float), -1.0, 1.0)


TET_4, TET_10, WEDGE_6, WEDGE_15, BRICK_8, BRICK_20 = (
    Tet4(), Tet10(), Wedge6(), Wedge15(), Brick8(), Brick20())

ELEMENTS = {e.name: e for e in (TET_4, TET_10, WEDGE_6, WEDGE_15, BRICK_8, BRICK_20)}

# node-count -> class, for face -> surface element lookup
SURFACE_TYPE_BY_NNODES = {3: "SURF_TRI_3", 4: "SURF_QUAD_4", 6: "SURF_TRI_6", 8: "SURF_QUAD_8"}


# ---------------------------------------------------------------------------
# surface (face) elements -- only what is needed for surface normals
# ---------------------------------------------------------------------------
def _tri3(r, s):
    return np.array([r, s, 1.0 - r - s])


def _tri6(r, s):
    u = 1.0 - r - s
    return np.array([(2 * r - 1) * r, (2 * s - 1) * s, (2 * u - 1) * u,
                     4.0 * r * s, 4.0 * u * s, 4.0 * u * r])


def _quad4(r, s):
    return np.array([(1 - r) * (1 - s), r * (1 - s), r * s, (1 - r) * s])


def _quad8(r, s):
    rr, rs, ss = r * r, r * s, s * s
    rrs, rss = rr * s, r * ss
    return np.array([0.25 * (rr - rrs + rs - rss + ss - 1.0),
                     0.25 * (rr - rrs - rs + rss + ss - 1.0),
                     0.25 * (rr + rrs + rs + rss + ss - 1.0),
                     0.25 * (rr + rrs - rs - rss + ss - 1.0),
                     0.5 * (-rr + rrs - s + 1.0),
                     0.5 * (r - rss - ss + 1.0),
                     0.5 * (-rr - rrs + s + 1.0),
                     0.5 * (-r + rss - ss + 1.0)])


# natural coordinates of the nodes of each surface element type
_SURF_NODE_RS = {
    3: (_tri3, np.array([[1, 0], [0, 1], [0, 0]], dtype=float)),
    4: (_quad4, np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)),
    6: (_tri6, np.array([[1, 0], [0, 1], [0, 0], [0.5, 0.5], [0, 0.5], [0.5, 0]])),
    8: (_quad8, np.array([[-1, -1], [1, -1], [1, 1], [-1, 1],
                          [0, -1], [1, 0], [0, 1], [-1, 0]], dtype=float)),
}


def surface_node_tangents(num_nodes, local_index, coords):
    """Tangent vectors (dx/dr, dx/ds) of a surface element at one of its nodes.

    ``coords`` is the (num_nodes, 3) array of the element's node coordinates in
    face order.  The (unnormalized) surface normal is ``cross(J0, J1)``.
    """
    shape_fn, node_rs = _SURF_NODE_RS[num_nodes]
    r, s = node_rs[local_index]
    d = np.empty((num_nodes, 2))
    d[:, 0] = shape_fn(complex(r, _H), s).imag / _H
    d[:, 1] = shape_fn(r, complex(s, _H)).imag / _H
    return d[:, 0] @ coords, d[:, 1] @ coords

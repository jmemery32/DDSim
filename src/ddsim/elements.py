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

The shape functions are implemented once, in :mod:`_kernels` (numba-compiled).
Their derivatives are obtained by complex-step differentiation: the shape
functions are polynomials, so this is exact to machine precision and removes
~300 lines of hand-derived derivative tables (and their chance of typos).
"""
import numpy as np

from . import _kernels as _K

_H = 1e-30  # complex-step size; exact for polynomials


# ---------------------------------------------------------------------------
# solid elements -- thin wrappers over the numba kernels in _kernels.py
# ---------------------------------------------------------------------------
class SolidElement:
    """Reference solid element (shape functions live in :mod:`_kernels`)."""

    name = ""
    code = -1
    num_nodes = 0
    center = (0.0, 0.0, 0.0)
    faces = ()          # tuple of node-index tuples (corner nodes first, then mid-side)

    def shape(self, nc):
        """Shape function values at natural coordinates ``nc`` -> ndarray(num_nodes)."""
        out = np.empty(_K.MAX_NODES)
        _K.shape_only(self.code, float(nc[0]), float(nc[1]), float(nc[2]), out,
                      np.empty(_K.MAX_NODES, dtype=np.complex128))
        return out[:self.num_nodes].copy()

    def dshape(self, nc):
        """Shape derivatives d N_k / d nc_j -> ndarray(num_nodes, 3)."""
        n = np.empty(_K.MAX_NODES)
        d = np.empty((_K.MAX_NODES, 3))
        _K.shape_dshape(self.code, float(nc[0]), float(nc[1]), float(nc[2]), n, d,
                        np.empty(_K.MAX_NODES, dtype=np.complex128))
        return d[:self.num_nodes].copy()

    def point_inside(self, nc, tol):
        return bool(_K.point_inside(self.code, np.asarray(nc, dtype=float), float(tol)))

    def nearest_point(self, nc):
        """Nearest point of the reference element to ``nc`` (in natural coordinates).

        For the simplex-type elements this is a true Euclidean projection; the 2007
        code moved the coordinates by ``-u/3`` with ``u = 1 - sum < 0`` (i.e. the wrong
        way) and, for the 10-node tet, mixed up r and s.  It only affects the
        reported distance for points that are (slightly) outside an element.
        """
        out = np.empty(3)
        _K.nearest_point(self.code, np.asarray(nc, dtype=float), out)
        return out


class Tet4(SolidElement):
    name, code, num_nodes = "TET_4", _K.TET_4, 4
    center = (0.25, 0.25, 0.25)
    faces = ((0, 1, 2), (1, 3, 2), (2, 3, 0), (0, 3, 1))


class Tet10(SolidElement):
    name, code, num_nodes = "TET_10", _K.TET_10, 10
    center = (0.25, 0.25, 0.25)
    faces = ((0, 1, 2, 4, 5, 6), (1, 3, 2, 8, 7, 5),
             (2, 3, 0, 7, 9, 6), (0, 3, 1, 9, 8, 4))


class Wedge6(SolidElement):
    name, code, num_nodes = "WEDGE_6", _K.WEDGE_6, 6
    center = (1.0 / 3.0, 1.0 / 3.0, 0.5)
    faces = ((0, 3, 4, 1), (2, 5, 3, 0), (1, 4, 5, 2), (0, 1, 2), (3, 5, 4))


class Wedge15(SolidElement):
    name, code, num_nodes = "WEDGE_15", _K.WEDGE_15, 15
    center = (1.0 / 3.0, 1.0 / 3.0, 0.0)
    faces = ((0, 3, 4, 1, 12, 9, 13, 6), (2, 5, 3, 0, 14, 11, 12, 8),
             (1, 4, 5, 2, 13, 10, 14, 7), (0, 1, 2, 6, 7, 8), (3, 5, 4, 11, 10, 9))


class Brick8(SolidElement):
    name, code, num_nodes = "BRICK_8", _K.BRICK_8, 8
    center = (0.5, 0.5, 0.5)
    faces = ((0, 1, 3, 2), (4, 6, 7, 5), (0, 4, 5, 1),
             (2, 3, 7, 6), (0, 2, 6, 4), (1, 5, 7, 3))


class Brick20(SolidElement):
    name, code, num_nodes = "BRICK_20", _K.BRICK_20, 20
    center = (0.0, 0.0, 0.0)
    faces = ((0, 1, 3, 2, 8, 13, 9, 12), (4, 6, 7, 5, 14, 11, 15, 10),
             (0, 4, 5, 1, 16, 10, 17, 8), (2, 3, 7, 6, 9, 19, 11, 18),
             (0, 2, 6, 4, 12, 18, 14, 16), (1, 5, 7, 3, 17, 15, 19, 13))


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

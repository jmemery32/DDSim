"""Geometry utilities for DDSim: ellipse (crack front) vs. surface-mesh crossings.

Pure-Python/numpy replacement for the compiled ``GeomUtils.pyd`` (built from
``GeomUtils.cpp`` / ``MeshWRangeTree.cpp``).  DDSim uses three entry points::

    cmesh  = GeomUtils.BuildSurfMeshCObject(facets)           # once per model
    thetas = GeomUtils.EllipseCMeshIntersections(cmesh, center, rot, a, b)
    length = GeomUtils.EllipseArcLength(a, b, theta_1, theta_2)

Ellipse convention (unchanged from 2007)
----------------------------------------
``rot`` is a flat 9-sequence (a row-major 3x3 matrix).  Its first row ``e1`` is
the direction of the semi-axis ``a`` and its second row ``e2`` that of ``b``::

    P(theta) = center + r(theta) * (cos(theta) e1 + sin(theta) e2)
    r(theta) = a b / sqrt(b^2 cos^2 + a^2 sin^2)             (polar form)

``theta`` is the *polar* angle (``phi``/"eccentric" angle: ``tan(psi) = (a/b)
tan(theta)``); the returned crossings are in ``[0, 2 pi)``.

How this differs from the C++
-----------------------------
* The ellipse/plane crossing is computed in closed form (the C++ used a
  3-point bracket + bisection that could miss thin intersections).
* Facets that share an edge report a crossing on that edge once, not twice.
* The tangency filter keeps its bookkeeping aligned (the C++ read ``ThCross``
  entries that were stale after removals) and prints nothing.
* Bounding-box candidate search pads by a small epsilon (the C++ padded every
  facet box by its full diagonal) -- the result set is unchanged.
* ``bisections`` / ``fraction`` are accepted for API compatibility and ignored.
"""
import math

import numpy as np

_TWO_PI = 2.0 * math.pi
_GL_X, _GL_W = np.polynomial.legendre.leggauss(9)


def _xyz(v):
    return np.array([v[0], v[1], v[2]], dtype=float)


# ---------------------------------------------------------------------------
# surface mesh
# ---------------------------------------------------------------------------
def _triangulate(corners):
    """Split a facet into triangles exactly as MeshWRangeTree did (quads along the
    shorter diagonal); facets with more corners are fanned around the centroid."""
    n = len(corners)
    if n == 3:
        return [corners]
    if n == 4:
        c = corners
        if np.linalg.norm(c[2] - c[0]) < np.linalg.norm(c[3] - c[1]):
            return [[c[0], c[1], c[2]], [c[0], c[2], c[3]]]
        return [[c[1], c[2], c[3]], [c[3], c[0], c[1]]]
    center = np.mean(corners, axis=0)
    return [[corners[j], corners[(j + 1) % n], center] for j in range(n)]


class SurfaceMesh:
    """A triangulated surface with per-triangle planes and bounding boxes."""

    def __init__(self, facets):
        tris = []
        for facet in facets:
            tris.extend(_triangulate([_xyz(v) for v in facet]))
        self.tri = np.array(tris, dtype=float).reshape(-1, 3, 3)
        e1 = self.tri[:, 1] - self.tri[:, 0]
        e2 = self.tri[:, 2] - self.tri[:, 0]
        n = np.cross(e1, e2)
        length = np.linalg.norm(n, axis=1)
        keep = length > 0.0                       # drop degenerate triangles
        self.tri, n, length = self.tri[keep], n[keep], length[keep]
        self.normal = n / length[:, None]
        self.offset = -np.einsum("ij,ij->i", self.tri[:, 1], self.normal)
        self.lo = self.tri.min(axis=1)
        self.hi = self.tri.max(axis=1)

    def __len__(self):
        return len(self.tri)

    def RangeQuery(self, lo, hi):
        """Indices of the triangles whose bounding boxes overlap [lo, hi]."""
        hit = np.all((self.lo <= hi) & (self.hi >= lo), axis=1)
        return np.nonzero(hit)[0]


def BuildSurfMeshCObject(facets):
    """``facets``: sequence of facets, each a sequence of 3 or 4 corner points."""
    return SurfaceMesh(facets)


# ---------------------------------------------------------------------------
# point in facet / arc length
# ---------------------------------------------------------------------------
def PointInFacet(point, facet, tol=0.0):
    """True if ``point`` (assumed to lie in the facet's plane) is inside the facet.

    Angle-sum test: the signed angles the polygon edges subtend at the point add
    up to 2 pi inside, pi on an edge and 0 outside; the threshold is pi/2.
    """
    p = _xyz(point)
    f = np.array([_xyz(v) for v in facet]) - p
    n = len(f)
    pt = np.zeros(3)
    if tol:
        ave = f.mean(axis=0)
        mag = np.linalg.norm(ave)
        if mag > 0.0:
            pt = ave / mag * tol
    a = f - pt
    b = np.roll(a, -1, axis=0)
    ma, mb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
    if np.any(ma == 0.0):                        # the point is a vertex
        return True
    # A point on an edge makes arccos() see (+-1 +- rounding) and its cross-product
    # sign is noise, which can cancel the sum.  Decide those explicitly.
    seg = b - a
    seg2 = np.einsum("ij,ij->i", seg, seg)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.clip(-np.einsum("ij,ij->i", a, seg) / seg2, 0.0, 1.0)
    closest = np.linalg.norm(a + t[:, None] * seg, axis=1)
    if np.any(closest <= 1e-10 * np.sqrt(seg2)):
        return True
    normal = np.cross(f[1] - f[0], f[2] - f[1])
    cosang = np.clip(np.einsum("ij,ij->i", a, b) / (ma * mb), -1.0, 1.0)
    sign = np.sign(np.einsum("j,ij->i", normal, np.cross(a, b)))
    return bool(abs(np.sum(sign * np.arccos(cosang))) > math.pi / 2.0)


def EllipseArcLength(a, b, th1, th2, dth=None):
    """Arc length of the ellipse between polar angles ``th1`` and ``th2``.

    Composite 9-point Gauss-Legendre on ``ds/dtheta = sqrt(r^2 + r'^2)`` (the
    2007 integrand).  ``dth`` is the largest sub-interval; the 2007 default of
    pi/2 leaves ~0.1% error for a 3:1 ellipse and much more for slender ones, so
    by default the step shrinks with the aspect ratio ``min(a,b)/max(a,b)``.
    """
    if dth is None:
        dth = min(0.2, max(0.005, 0.2 * min(a, b) / max(a, b)))
    span = th2 - th1
    numstep = max(1, int(math.ceil(abs(span) / dth)))
    step = span / numstep
    edges = th1 + step * np.arange(numstep + 1)
    lo, hi = edges[:-1, None], edges[1:, None]
    th = 0.5 * (hi + lo) + 0.5 * (hi - lo) * _GL_X
    a2, b2 = a * a, b * b
    s2, c2 = np.sin(th) ** 2, np.cos(th) ** 2
    q = a2 * s2 + b2 * c2
    f = a2 * b2 * (a2 - b2) ** 2 * s2 * c2 / q ** 3 + a2 * b2 / q
    return abs(float(np.sum(_GL_W * np.sqrt(f)))) * abs(step) / 2.0


# ---------------------------------------------------------------------------
# ellipse / plane
# ---------------------------------------------------------------------------
def _frame(rot):
    r = np.asarray(rot, dtype=float).reshape(-1)
    if r.size != 9:
        raise TypeError("argument does not match expected size")
    return r[0:3], r[3:6]


def _polar_angle(a, b, psi):
    """Polar angle theta in [0, 2 pi) of the point (a cos psi, b sin psi)."""
    return math.atan2(b * math.sin(psi), a * math.cos(psi)) % _TWO_PI


def _plane_crossings(p, q, d0, a, b):
    """Polar angles where ``d0 + p*x + q*y = 0`` meets the ellipse x^2/a^2 + y^2/b^2 = 1.

    Returns ``(theta_up, theta_down)`` -- the crossing where the plane distance
    increases with theta, then the one where it decreases -- or ``None``.
    """
    amp = math.hypot(p * a, q * b)
    if amp == 0.0:
        return None                               # plane parallel to the ellipse plane
    c = -d0 / amp
    if abs(c) > 1.0:
        if abs(c) > 1.0 + 1e-12:
            return None
        c = math.copysign(1.0, c)
    psi0 = math.atan2(q * b, p * a)
    alpha = math.acos(c)
    return _polar_angle(a, b, psi0 - alpha), _polar_angle(a, b, psi0 + alpha)


def EllipsePlaneIntersection(normal, plane_d, center, rot, a, b, bisections=14):
    """Intersection of the ellipse with the plane ``n . x + d = 0`` (n is normalised,
    ``d`` is used as given).  Returns the two polar angles ``(theta_1, theta_2)``,
    or ``(-1, -1)`` if the ellipse does not meet the plane."""
    n = _xyz(normal)
    n = n / np.linalg.norm(n)
    e1, e2 = _frame(rot)
    d0 = float(n @ _xyz(center)) + plane_d
    res = _plane_crossings(float(n @ e1), float(n @ e2), d0, a, b)
    return (-1, -1) if res is None else res


# ---------------------------------------------------------------------------
# ellipse / mesh
# ---------------------------------------------------------------------------
def _ellipse_point(center, e1, e2, a, b, theta):
    r = a * b / math.sqrt(b * b * math.cos(theta) ** 2 + a * a * math.sin(theta) ** 2)
    return center + r * (math.cos(theta) * e1 + math.sin(theta) * e2)


def _crossings(mesh, tri_ids, center, e1, e2, a, b):
    """[(theta, j, triangle)] for all crossings of the ellipse with the triangles."""
    hits = []
    for t in tri_ids:
        n, d = mesh.normal[t], mesh.offset[t]
        res = _plane_crossings(float(n @ e1), float(n @ e2), float(n @ center) + d, a, b)
        if res is None:
            continue
        pair = sorted(res)                        # j = 0: smaller theta, 1: larger
        for j, theta in enumerate(pair):
            if PointInFacet(_ellipse_point(center, e1, e2, a, b, theta), mesh.tri[t]):
                hits.append((theta, j, int(t)))
    return hits


def _remove_tangencies(hits, a, b):
    """The 2007 sliver filter.

    Crossings closer together (along the ellipse) than 10% of its perimeter --
    or, seen the other way round, farther apart than 90% -- are treated as a
    tangent touch rather than a real excursion: an adjacent (0, 1) pair is dropped
    completely, otherwise only the first of the two.
    """
    h = (a - b) ** 2 / (a + b) ** 2
    p = math.pi * (a + b) * (1.0 + 3.0 * h / (10.0 + math.sqrt(4.0 - 3.0 * h)))   # Ramanujan

    def is_sliver(t0, t1):
        arc = EllipseArcLength(a, b, t0, t1)
        return abs(1.0 - arc / p) <= 0.1 or abs(1.0 - arc / p) >= 0.9

    th = [hit[0] for hit in hits]
    cr = [hit[1] for hit in hits]
    alive = [True] * len(th)

    i = 0
    while i < len(th) - 1:
        if is_sliver(th[i], th[i + 1]):
            if cr[i] + cr[i + 1] == 1:
                alive[i] = alive[i + 1] = False
                i += 1
            else:
                alive[i] = False
        i += 1
    th = [t for t, ok in zip(th, alive) if ok]
    cr = [c for c, ok in zip(cr, alive) if ok]

    if len(th) > 1 and is_sliver(th[0], th[-1]):
        if cr[0] + cr[-1] == 1:
            th = th[1:-1]
        else:
            th = th[1:]
    return th


def _ellipse_mesh_intersections(mesh, center, rot, a, b):
    center = _xyz(center)
    e1, e2 = _frame(rot)

    # the ellipse lies inside the bounding box of its enclosing rectangle
    corners = np.array([center + sa * a * e1 + sb * b * e2 for sa in (1, -1) for sb in (1, -1)])
    lo, hi = corners.min(axis=0), corners.max(axis=0)
    eps = 1e-9 * max(a, b)
    ids = mesh.RangeQuery(lo - eps, hi + eps)

    hits = _crossings(mesh, ids, center, e1, e2, a, b)
    hits.sort(key=lambda h: h[0])                 # stable, like the C++ bubble sort

    # the same point on an edge shared by two triangles is one crossing
    merged = []
    for h in hits:
        if merged and merged[-1][2] != h[2] and abs(h[0] - merged[-1][0]) < 1e-9:
            continue
        merged.append(h)
    if (len(merged) > 1 and merged[0][2] != merged[-1][2]
            and abs(merged[0][0] + _TWO_PI - merged[-1][0]) < 1e-9):
        merged.pop()                              # the pair straddles theta = 0

    return _remove_tangencies(merged, a, b)


def EllipseCMeshIntersections(mesh, center, rot, a, b, bisections=8, fraction=0.1):
    """Polar angles (sorted, in [0, 2 pi)) where the ellipse crosses the surface mesh."""
    return _ellipse_mesh_intersections(mesh, center, rot, a, b)


def EllipseMeshIntersections(facets, center, rot, a, b, bisections=8, fraction=0.1):
    """Same as :func:`EllipseCMeshIntersections` for a plain list of facets."""
    return _ellipse_mesh_intersections(SurfaceMesh(facets), center, rot, a, b)

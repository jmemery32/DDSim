import math

import numpy as np
import pytest
from scipy import integrate
from scipy.special import ellipe

from ddsim import GeomUtils as G
from ddsim import Vec3D

I3 = [1, 0, 0, 0, 1, 0, 0, 0, 1]


def cube_facets(h):
    """Six quad facets of the cube [-h, h]^3 (corner Vec3Ds, like DamMo builds)."""
    v = lambda x, y, z: Vec3D.Vec3D(x, y, z)
    return [
        [v(-h, -h, h), v(h, -h, h), v(h, h, h), v(-h, h, h)],
        [v(-h, -h, -h), v(-h, h, -h), v(h, h, -h), v(h, -h, -h)],
        [v(h, -h, -h), v(h, h, -h), v(h, h, h), v(h, -h, h)],
        [v(-h, -h, -h), v(-h, -h, h), v(-h, h, h), v(-h, h, -h)],
        [v(-h, h, -h), v(-h, h, h), v(h, h, h), v(h, h, -h)],
        [v(-h, -h, -h), v(h, -h, -h), v(h, -h, h), v(-h, -h, h)],
    ]


def rot_z(angle):
    c, s = math.cos(angle), math.sin(angle)
    return [c, s, 0, -s, c, 0, 0, 0, 1]


# ---------------------------------------------------------------------------
@pytest.mark.parametrize("a,b", [(3.0, 1.0), (1.0, 3.0), (10.0, 1.0), (1.0, 1.0), (2.5, 2.4)])
def test_arc_length_full_ellipse_matches_elliptic_integral(a, b):
    exact = 4.0 * max(a, b) * ellipe(1.0 - (min(a, b) / max(a, b)) ** 2)
    assert G.EllipseArcLength(a, b, 0.0, 2 * math.pi) == pytest.approx(exact, rel=1e-9)


def test_2007_default_step_was_inaccurate_for_slender_ellipses():
    # documents why the default step now depends on the aspect ratio
    exact = 4.0 * 3.0 * ellipe(1.0 - 1.0 / 9.0)
    coarse = G.EllipseArcLength(3.0, 1.0, 0.0, 2 * math.pi, 1.570796)
    assert abs(coarse / exact - 1.0) > 1e-3


@pytest.mark.parametrize("a,b", [(3.0, 1.0), (1.0, 2.5), (1.0, 1.0)])
def test_arc_length_matches_independent_quadrature(a, b):
    def ds(t):
        r = lambda t: a * b / math.sqrt(b * b * math.cos(t) ** 2 + a * a * math.sin(t) ** 2)
        x = lambda t: r(t) * math.cos(t)
        y = lambda t: r(t) * math.sin(t)
        h = 1e-6
        return math.hypot((x(t + h) - x(t - h)) / (2 * h), (y(t + h) - y(t - h)) / (2 * h))

    for t1, t2 in [(0.2, 1.0), (0.0, 3.0), (1.0, 5.5)]:
        exact, _ = integrate.quad(ds, t1, t2)
        assert G.EllipseArcLength(a, b, t1, t2) == pytest.approx(exact, rel=1e-6)


def test_arc_length_is_additive_and_symmetric():
    whole = G.EllipseArcLength(2, 1, 0.3, 4.0)
    assert whole == pytest.approx(G.EllipseArcLength(2, 1, 0.3, 2.0) + G.EllipseArcLength(2, 1, 2.0, 4.0))
    assert G.EllipseArcLength(2, 1, 4.0, 0.3) == pytest.approx(whole)


# ---------------------------------------------------------------------------
@pytest.mark.parametrize("a,b", [(2.0, 1.0), (1.0, 2.0), (1.5, 1.5)])
def test_plane_intersection_is_on_both_ellipse_and_plane(a, b):
    rng = np.random.default_rng(5)
    rot = rot_z(0.4)
    center = np.array([0.2, -0.1, 0.3])
    e1, e2 = np.array(rot[0:3]), np.array(rot[3:6])
    found = 0
    for _ in range(50):
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        d = rng.uniform(-1.5, 1.5)
        res = G.EllipsePlaneIntersection(n, d, center, rot, a, b)
        if res == (-1, -1):
            continue
        found += 1
        for th in res:
            r = a * b / math.sqrt(b * b * math.cos(th) ** 2 + a * a * math.sin(th) ** 2)
            pt = center + r * (math.cos(th) * e1 + math.sin(th) * e2)
            assert n @ pt + d == pytest.approx(0.0, abs=1e-10)
    assert found > 10


def test_plane_intersection_none_when_plane_misses():
    assert G.EllipsePlaneIntersection([0, 0, 1], -100.0, [0, 0, 0], I3, 2, 1) == (-1, -1)
    assert G.EllipsePlaneIntersection([0, 0, 1], 0.5, [0, 0, 0], I3, 2, 1) == (-1, -1)  # parallel


def test_plane_intersection_up_then_down():
    n, d = [1.0, 0.0, 0.0], -1.0                 # plane x = 1, ellipse x semi-axis 2
    up, down = G.EllipsePlaneIntersection(n, d, [0, 0, 0], I3, 2, 1)
    dist = lambda th: (2 * math.cos(th) * (1 / math.sqrt(1)) if False else 0)
    # D(psi) = 2 cos(psi) - 1: falls through zero at psi = +pi/3 and rises at -pi/3
    th_of = lambda psi: math.atan2(math.sin(psi), 2 * math.cos(psi)) % (2 * math.pi)
    assert down == pytest.approx(th_of(math.pi / 3))
    assert up == pytest.approx(th_of(-math.pi / 3))


# ---------------------------------------------------------------------------
def test_point_in_facet():
    tri = [(0, 0, 0), (2, 0, 0), (0, 2, 0)]
    assert G.PointInFacet((0.5, 0.5, 0), tri)
    assert not G.PointInFacet((1.5, 1.5, 0), tri)
    assert not G.PointInFacet((-0.1, 0.5, 0), tri)
    assert G.PointInFacet((1.0, 0.0, 0), tri)                    # on an edge
    assert G.PointInFacet((0, 0, 0), tri)                        # on a vertex
    quad = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
    assert G.PointInFacet((0.5, 0.5, 0), quad)
    assert not G.PointInFacet((1.2, 0.5, 0), quad)
    reversed_tri = tri[::-1]
    assert G.PointInFacet((0.5, 0.5, 0), reversed_tri)           # orientation independent


def test_quads_are_split_into_triangles():
    mesh = G.BuildSurfMeshCObject(cube_facets(1.0))
    assert len(mesh) == 12
    assert mesh.tri.shape == (12, 3, 3)


# ---------------------------------------------------------------------------
def crossing_angle(a, b, x_plane):
    """theta of the crossing of x^2/a^2 + y^2/b^2 = 1 with x = x_plane, y > 0."""
    y = b * math.sqrt(1 - (x_plane / a) ** 2)
    return math.atan2(y, x_plane)


def test_ellipse_crossing_two_faces_gives_four_angles():
    # ellipse a=2, b=1 in the z=0 plane, cube half-width 1.5: it pokes out through x = +-1.5
    mesh = G.BuildSurfMeshCObject(cube_facets(1.5))
    th = G.EllipseCMeshIntersections(mesh, [0, 0, 0], I3, 2.0, 1.0)
    t0 = crossing_angle(2.0, 1.0, 1.5)
    expected = sorted([t0, math.pi - t0, math.pi + t0, 2 * math.pi - t0])
    assert th == pytest.approx(expected, abs=1e-9)


def test_ellipse_inside_mesh_has_no_crossings():
    mesh = G.BuildSurfMeshCObject(cube_facets(5.0))
    assert G.EllipseCMeshIntersections(mesh, [0, 0, 0], I3, 2.0, 1.0) == []


def test_swapped_axes_and_rotation():
    # a < b: the long axis is e2; rotate the frame 90 degrees so the long axis is x again
    mesh = G.BuildSurfMeshCObject(cube_facets(1.5))
    th = G.EllipseCMeshIntersections(mesh, [0, 0, 0], rot_z(math.pi / 2), 1.0, 2.0)
    assert len(th) == 4
    # rotating the frame by 90 degrees moves every polar angle back by 90 degrees
    t0 = crossing_angle(2.0, 1.0, 1.5)
    ref = [t0, math.pi - t0, math.pi + t0, 2 * math.pi - t0]
    expected = sorted(((t - math.pi / 2) % (2 * math.pi)) for t in ref)
    assert th == pytest.approx(expected, abs=1e-9)


def test_crossing_on_a_shared_edge_is_reported_once():
    # a unit-ish square in z = 0, split along the diagonal x + y = 0.  The ellipse lies
    # in the vertical plane y = -x and crosses z = 0 exactly ON that diagonal, i.e. on
    # the edge shared by the two triangles, at theta = 0 and theta = pi.
    square = [[Vec3D.Vec3D(-1, -1, 0), Vec3D.Vec3D(1, -1, 0), Vec3D.Vec3D(1, 1, 0), Vec3D.Vec3D(-1, 1, 0)]]
    mesh = G.BuildSurfMeshCObject(square)
    assert len(mesh) == 2
    s = 1 / math.sqrt(2.0)
    rot = [s, -s, 0, 0, 0, 1, s, s, 0]
    th = G.EllipseCMeshIntersections(mesh, [0, 0, 0], rot, 0.5, 0.3)
    assert len(th) == 2
    wrapped = sorted(t if t < math.pi + 0.5 else t - 2 * math.pi for t in th)
    assert wrapped == pytest.approx([0.0, math.pi], abs=1e-9)


def test_short_excursions_near_a_corner_are_filtered():
    # A circle just larger than the cube's half-diagonal leaves through each vertical edge
    # region and comes straight back: pairs of crossings a few % of the perimeter
    # apart.  The 2007 sliver filter treats those as noise.
    r = 1.5 * math.sqrt(2.0)
    mesh = G.BuildSurfMeshCObject(cube_facets(1.5))
    assert G.EllipseCMeshIntersections(mesh, [0, 0, 0], I3, r + 1e-3, r + 1e-3) == []


def test_tangent_touch_is_filtered():
    # circle radius 1.5 just touches the four faces of the cube at their centres
    mesh = G.BuildSurfMeshCObject(cube_facets(1.5))
    th = G.EllipseCMeshIntersections(mesh, [0, 0, 0], I3, 1.5, 1.5)
    assert th == []


def test_offset_center_and_tilted_plane_matches_brute_force():
    """3D case: compare with a dense sampling of the ellipse for inside/outside flips."""
    rng = np.random.default_rng(2)
    h = 1.5
    mesh = G.BuildSurfMeshCObject(cube_facets(h))
    for _ in range(20):
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(q) < 0:
            q[:, 0] = -q[:, 0]
        rot = list(q.T.reshape(-1))
        e1, e2 = q[:, 0], q[:, 1]
        center = rng.uniform(-0.5, 0.5, 3)
        a, b = rng.uniform(0.6, 2.4), rng.uniform(0.6, 2.4)
        th = G.EllipseCMeshIntersections(mesh, center, rot, a, b)

        theta = np.linspace(0, 2 * np.pi, 20001)
        r = a * b / np.sqrt(b * b * np.cos(theta) ** 2 + a * a * np.sin(theta) ** 2)
        pts = center + (r * np.cos(theta))[:, None] * e1 + (r * np.sin(theta))[:, None] * e2
        inside = np.all(np.abs(pts) <= h, axis=1)
        flips = np.nonzero(inside[:-1] != inside[1:])[0]
        brute = sorted(theta[i] for i in flips)
        # the 2007 filter may legitimately drop slivers; when nothing was dropped they agree
        if len(th) == len(brute):
            assert th == pytest.approx(brute, abs=2 * (theta[1] - theta[0]))
        else:
            assert len(th) < len(brute)


# ---------------------------------------------------------------------------
# compiled kernels vs the readable numpy reference
# ---------------------------------------------------------------------------
def test_point_in_triangle_kernel_matches_reference():
    from ddsim import _geom_kernels as GK
    rng = np.random.default_rng(9)
    for _ in range(300):
        tri = rng.normal(size=(3, 3))
        # points in the triangle's plane: inside, outside, on edges, on vertices
        w = rng.dirichlet([1, 1, 1])
        candidates = [w @ tri, tri[0], tri[1], 0.5 * (tri[0] + tri[1]), 0.5 * (tri[1] + tri[2]),
                      tri[0] + 3 * (tri[1] - tri[0]), (w * [2.5, -1.0, -0.5]) @ tri]
        for p in candidates:
            assert GK.point_in_triangle(np.array(p), tri) == G.PointInFacet(p, tri)


def test_crossings_kernel_matches_reference():
    rng = np.random.default_rng(12)
    total = 0
    for trial in range(40):
        mesh = G.BuildSurfMeshCObject(cube_facets(rng.uniform(0.8, 1.8)))
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        e1, e2 = q[:, 0], q[:, 1]
        center = rng.uniform(-0.4, 0.4, 3)
        a, b = rng.uniform(0.5, 2.6), rng.uniform(0.5, 2.6)
        ids = np.arange(len(mesh))
        fast = G._crossings(mesh, ids, center, e1, e2, a, b)
        ref = G._crossings_reference(mesh, ids, center, e1, e2, a, b)
        assert [(j, t) for _, j, t in fast] == [(j, t) for _, j, t in ref]
        assert [h[0] for h in fast] == pytest.approx([h[0] for h in ref], abs=1e-12)
        total += len(ref)
    assert total > 40                       # the comparison actually saw crossings

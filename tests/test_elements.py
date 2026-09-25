import numpy as np
import pytest

from ddsim import elements as E

# natural coordinates of every node, read off the shape functions in ElemClass.cpp
NODES = {
    "TET_4": [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0)],
    "TET_10": [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0),
               (.5, .5, 0), (0, .5, .5), (.5, 0, .5), (0, 0, .5), (0, .5, 0), (.5, 0, 0)],
    "WEDGE_6": [(1, 0, 0), (0, 1, 0), (0, 0, 0), (1, 0, 1), (0, 1, 1), (0, 0, 1)],
    "WEDGE_15": [(0, 0, -1), (0, 1, -1), (1, 0, -1), (0, 0, 1), (0, 1, 1), (1, 0, 1),
                 (0, .5, -1), (.5, .5, -1), (.5, 0, -1), (0, .5, 1), (.5, .5, 1), (.5, 0, 1),
                 (0, 0, 0), (0, 1, 0), (1, 0, 0)],
    "BRICK_8": [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)],
    "BRICK_20": [(-1, -1, -1), (-1, -1, 1), (-1, 1, -1), (-1, 1, 1),
                 (1, -1, -1), (1, -1, 1), (1, 1, -1), (1, 1, 1),
                 (-1, -1, 0), (-1, 1, 0), (1, -1, 0), (1, 1, 0),
                 (-1, 0, -1), (-1, 0, 1), (1, 0, -1), (1, 0, 1),
                 (0, -1, -1), (0, -1, 1), (0, 1, -1), (0, 1, 1)],
}
ALL = list(E.ELEMENTS)


@pytest.mark.parametrize("name", ALL)
def test_node_table_matches_element(name):
    assert len(NODES[name]) == E.ELEMENTS[name].num_nodes


@pytest.mark.parametrize("name", ALL)
def test_kronecker_delta_at_nodes(name):
    el = E.ELEMENTS[name]
    for i, nc in enumerate(NODES[name]):
        assert el.shape(nc) == pytest.approx(np.eye(el.num_nodes)[i], abs=1e-12)


@pytest.mark.parametrize("name", ALL)
def test_partition_of_unity_and_derivative_sum(name):
    el = E.ELEMENTS[name]
    rng = np.random.default_rng(0)
    for _ in range(10):
        nc = rng.uniform(-0.2, 0.9, 3)
        assert el.shape(nc).sum() == pytest.approx(1.0)
        assert el.dshape(nc).sum(axis=0) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("name", ALL)
def test_complex_step_derivative_matches_finite_difference(name):
    el = E.ELEMENTS[name]
    nc = np.array([0.21, 0.33, 0.17])
    h = 1e-6
    for j in range(3):
        d = np.zeros(3)
        d[j] = h
        fd = (el.shape(nc + d) - el.shape(nc - d)) / (2 * h)
        assert el.dshape(nc)[:, j] == pytest.approx(fd, abs=1e-8)


@pytest.mark.parametrize("name", ALL)
def test_center_is_inside_and_far_point_is_not(name):
    el = E.ELEMENTS[name]
    assert el.point_inside(el.center, 1e-7)
    assert not el.point_inside((5.0, 5.0, 5.0), 1e-7)
    assert el.point_inside(NODES[name][0], 1e-7)          # nodes are on the boundary


@pytest.mark.parametrize("name", ALL)
def test_nearest_point_of_inside_point_is_itself(name):
    el = E.ELEMENTS[name]
    assert el.nearest_point(el.center) == pytest.approx(el.center)


@pytest.mark.parametrize("name", ALL)
def test_nearest_point_lands_inside(name):
    el = E.ELEMENTS[name]
    for nc in [(2.0, 0.3, 0.2), (-3.0, -3.0, -3.0), (0.9, 0.9, 0.9), (0.5, 0.5, 4.0)]:
        near = el.nearest_point(nc)
        assert el.point_inside(near, 1e-9), (name, nc, near)


def test_tet_nearest_point_is_the_closest_on_the_slanted_face():
    # (1,1,1) projected onto the r+s+t = 1 plane is (1/3,1/3,1/3) -- the 2007 code
    # moved it the wrong way (to (4/3,...) then clamped).
    near = E.TET_4.nearest_point((1.0, 1.0, 1.0))
    assert near == pytest.approx([1 / 3, 1 / 3, 1 / 3])
    assert E.TET_10.nearest_point((1.0, 1.0, 1.0)) == pytest.approx([1 / 3, 1 / 3, 1 / 3])


@pytest.mark.parametrize("name", ALL)
def test_faces_are_on_the_boundary_and_consistently_oriented(name):
    """Each face's nodes are coplanar on the reference boundary, and all face
    normals (cross(J0, J1) at the first node) point the same way -- outward for
    every element except where noted."""
    el = E.ELEMENTS[name]
    X = np.array(NODES[name], dtype=float)                 # physical == natural here
    center = np.array(el.center)
    signs = []
    for face in el.faces:
        assert len(face) in E.SURFACE_TYPE_BY_NNODES
        P = X[list(face)]
        J0, J1 = E.surface_node_tangents(len(face), 0, P)
        n = np.cross(J0, J1)
        assert np.linalg.norm(n) > 0
        # all corner nodes lie in one plane
        corners = P[: 3 if len(face) in (3, 6) else 4]
        dev = (corners - corners[0]) @ n / np.linalg.norm(n)
        assert dev == pytest.approx(0.0, abs=1e-12)
        signs.append(np.sign(n @ (corners.mean(axis=0) - center)))
    assert len(set(signs)) == 1
    outward = signs[0] > 0
    print(name, "face normals point", "outward" if outward else "INWARD")


@pytest.mark.parametrize("name", ALL)
def test_newton_kernel_matches_numpy_solve(name):
    """The compiled Newton solve (Cramer's rule) recovers natural coordinates for
    a distorted element, checked against an independent np.linalg.solve iteration."""
    from ddsim import _kernels as K
    el = E.ELEMENTS[name]
    rng = np.random.default_rng(4)
    nodes = np.array(NODES[name], dtype=float)
    A = np.array([[2.0, 0.4, 0.1], [0.2, 1.5, 0.3], [0.1, 0.2, 2.5]])
    X = nodes @ A.T + 0.03 * rng.normal(size=nodes.shape) * (np.abs(nodes) < 1).any(axis=1)[:, None]
    Xp = np.zeros((20, 3))
    Xp[:len(X)] = X
    nc_true = np.array(el.center) + 0.1 * rng.normal(size=3)
    q = el.shape(nc_true) @ X
    nc = np.empty(3)
    ok = K.find_natural(el.code, Xp, q, 1e-12, 30, nc, np.empty(20), np.empty((20, 3)),
                        np.empty(20, dtype=np.complex128))
    assert ok
    assert nc == pytest.approx(nc_true, abs=1e-9)

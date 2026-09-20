import os

import numpy as np
import pytest

import MeshTools
import Vec3D
import elements as E
import mesh_io
from test_elements import NODES

EXAMPLE1 = os.path.join(os.path.dirname(__file__), "..", "LevelIIIv0.1", "models", "example1")

A = np.array([[2.0, 0.3, 0.0], [0.1, 1.5, 0.2], [0.0, 0.4, 3.0]])   # det > 0, affine map
B = np.array([1.0, 2.0, 3.0])


def single_element_mesh(name):
    """One element of the given type, geometry x = A nc + B, stress linear in x."""
    nc = np.array(NODES[name], dtype=float)
    xyz = nc @ A.T + B
    G = np.random.default_rng(7).normal(size=(6, 3))
    c0 = np.arange(1.0, 7.0)
    stress = c0 + xyz @ G.T
    data = mesh_io.MeshData(
        node_ids=np.arange(1, len(nc) + 1), coords=xyz,
        elem_ids=np.array([10]), elem_types=[name],
        connectivity=[np.arange(1, len(nc) + 1)], stress=stress)
    field = lambda x: c0 + G @ x
    return MeshTools.MeshTools(data=data), field, nc


def as_vec(x):
    return Vec3D.Vec3D(*x)


# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", list(E.ELEMENTS))
def test_linear_field_is_reproduced_exactly(name):
    model, field, _ = single_element_mesh(name)
    el = E.ELEMENTS[name]
    rng = np.random.default_rng(3)
    for _ in range(15):
        # random point strictly inside the reference element
        while True:
            nc = rng.uniform(-1, 1.0, 3)
            if name in ("BRICK_8", "TET_4", "TET_10", "WEDGE_6"):
                nc = np.abs(nc)
            if el.point_inside(nc, -1e-3):
                break
        x = A @ nc + B
        sig = model.GetPtStress(as_vec(x))
        got = np.array([sig.xx(), sig.yy(), sig.zz(), sig.xy(), sig.yz(), sig.zx()])
        assert got == pytest.approx(field(x), rel=1e-8, abs=1e-8)
        assert model.IsPointOutsideMesh(as_vec(x))[0] == -2


@pytest.mark.parametrize("name", list(E.ELEMENTS))
def test_far_point_is_outside_and_raises(name):
    model, _, _ = single_element_mesh(name)
    far = as_vec(np.array([1e3, 1e3, 1e3]))
    assert model.IsPointOutsideMesh(far) == (-3, 0.0)
    with pytest.raises(MeshTools.EmptySearchResult):
        model.GetPtStress(far)


@pytest.mark.parametrize("name", list(E.ELEMENTS))
def test_just_outside_point_reports_distance(name):
    model, _, nc = single_element_mesh(name)
    corner = A @ nc[0] + B
    centroid = (nc @ A.T + B).mean(axis=0)
    out = corner + 1e-3 * (corner - centroid) / np.linalg.norm(corner - centroid)
    status, dist = model.IsPointOutsideMesh(as_vec(out))
    assert status == -1
    assert 0.0 < dist < 5e-3


@pytest.mark.parametrize("name", list(E.ELEMENTS))
def test_surface_normals_point_outward(name):
    model, _, nc = single_element_mesh(name)
    xyz = nc @ A.T + B
    centroid = xyz.mean(axis=0)
    for nid in model.GetNodeList():
        if not model.IsSurfaceNode(nid):
            continue
        n = model.SurfaceNormal(nid)
        assert n.Magnitude() == pytest.approx(1.0)
        assert n * as_vec(xyz[nid - 1] - centroid) > 0.0, (name, nid)
        for seid in model.GetAdjacentSurfElems(nid):
            assert model.SurfaceNormal(nid, seid).Magnitude() == pytest.approx(1.0)


def brick_stack(n):
    """n Brick8 elements stacked in x -- n*4 + 4 nodes."""
    nodes, coords = {}, []
    for i in range(n + 1):
        for j in (0, 1):
            for k in (0, 1):
                nodes[(i, j, k)] = len(coords) + 1
                coords.append((float(i), float(j), float(k)))
    conn, ids = [], []
    for e in range(n):
        c = [nodes[(e + a, b, cc)] for a in (0, 1) for b in (0, 1) for cc in (0, 1)]
        # Brick8 order: r is the slowest index, then s, then t
        conn.append(np.array(c))
        ids.append(e + 1)
    data = mesh_io.MeshData(node_ids=np.arange(1, len(coords) + 1), coords=np.array(coords),
                            elem_ids=np.array(ids), elem_types=["BRICK_8"] * n,
                            connectivity=conn, stress=np.zeros((len(coords), 6)))
    return MeshTools.MeshTools(data=data)


def test_shared_face_is_not_on_the_surface():
    m = brick_stack(2)
    assert len(m._surf_nodes) == 10                    # 12 faces - the 2 shared ones
    assert all(m.IsSurfaceNode(n) for n in m.GetNodeList())
    assert m.GetAdjacentElems(5) == [1, 2]             # middle nodes touch both bricks
    assert m.GetElemInfo(2) == [5, 6, 7, 8, 9, 10, 11, 12]
    with pytest.raises(MeshTools.InvalidNodeId):
        m.GetAdjacentElems(999)
    with pytest.raises(MeshTools.InvalidElemId):
        m.GetElemInfo(999)


def test_point_on_shared_face_is_found_from_either_side():
    m = brick_stack(2)
    assert m.IsPointOutsideMesh(as_vec([1.0, 0.5, 0.5]))[0] == -2
    assert m.IsPointOutsideMesh(as_vec([1.0 + 1e-9, 0.5, 0.5]))[0] == -2
    assert m.IsPointOutsideMesh(as_vec([0.5, 0.5, 0.5]))[0] == -2      # exercises the cache
    assert m.IsPointOutsideMesh(as_vec([1.5, 0.5, 0.5]))[0] == -2


# ---------------------------------------------------------------------------
# the real example model (12 quadratic tets, uniform syy = 3.4)
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def example1():
    m = MeshTools.MeshTools(EXAMPLE1, "RDB")
    m.SetPointInsideTolerance(1.0e-7)
    return m


def test_example1_reads(example1):
    assert len(example1.GetNodeList()) == 35
    assert len(example1._elem_ids) == 12
    xyz, disp, sig = example1.GetNodeInfo(0)
    assert list(xyz) == [-10, 10, 10]
    assert list(disp) == [0, 0, 0]
    assert sig.yy() == pytest.approx(3.4)
    with pytest.raises(MeshTools.InvalidNodeId):
        example1.GetNodeInfo(12345)
    Max, Min = example1.GetMaxDimension()
    assert list(Max) == [10, 10, 10] and list(Min) == [-10, -10, -10]


def test_example1_uniform_stress_everywhere(example1):
    rng = np.random.default_rng(11)
    for _ in range(200):
        p = rng.uniform(-9.99, 9.99, 3)
        s = example1.GetPtStress(as_vec(p))
        assert (s.xx(), s.yy(), s.zz(), s.xy(), s.yz(), s.zx()) == pytest.approx(
            (0, 3.4, 0, 0, 0, 0), abs=1e-6)


def test_example1_surface(example1):
    surface = [n for n in example1.GetNodeList() if example1.IsSurfaceNode(n)]
    assert len(surface) == 26                          # 8 corners + 12 edge mids + 6 face centers
    assert not example1.IsSurfaceNode(8)               # the center node
    assert len(example1._surf_nodes) == 12 * 1         # 12 boundary triangles (6 faces x 2)

    n = example1.SurfaceNormal(10)                     # face-center node (0, 0, 10)
    assert (n.x(), n.y(), n.z()) == pytest.approx((0, 0, 1), abs=1e-9)
    n = example1.SurfaceNormal(12)                     # (-10, 0, 0)
    assert (n.x(), n.y(), n.z()) == pytest.approx((-1, 0, 0), abs=1e-9)
    corner = example1.SurfaceNormal(0)                 # (-10, 10, 10): mixes three faces
    assert corner * as_vec([-1, 1, 1]) / 3 ** 0.5 > 0.5


def test_example1_outside(example1):
    assert example1.IsPointOutsideMesh(as_vec([0, 0, 0]))[0] == -2
    status, dist = example1.IsPointOutsideMesh(as_vec([0, 0, 10.5]))
    assert status == -1 and dist == pytest.approx(0.5, abs=1e-3)
    assert example1.IsPointOutsideMesh(as_vec([0, 0, 500]))[0] == -3

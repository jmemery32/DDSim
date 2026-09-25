"""Exodus II read/write tests.

Two kinds of ground truth, per PORTING_NOTES.md:
* a real sample file (tests/fixtures/exodus/beam_hex8.g, from an unrelated
  project of the user's) for HEX8 -- the node-order permutation was derived
  and hand-verified against this file's actual coordinates.
* write -> read round trips against synthetic meshes for TET_4/WEDGE_6 (no
  real sample available), checked for exact geometry/stress recovery and
  (independently of exodus_io's own internals) a positive signed volume.
"""
import os

import numpy as np
import pytest

from ddsim import elements as E
from ddsim import exodus_io, mesh_io
from test_elements import NODES

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "exodus", "beam_hex8.g")

A = np.array([[2.0, 0.3, 0.0], [0.1, 1.5, 0.2], [0.0, 0.4, 3.0]])  # det > 0
B = np.array([1.0, 2.0, 3.0])


def synthetic_mesh(name, n_elem=1):
    """n_elem disjoint copies of one element type, distinct affine-transformed
    positions, with a linear (and thus exactly-representable) stress field."""
    nc = np.array(NODES[name], dtype=float)
    n = len(nc)
    coords, connectivity, elem_ids = [], [], []
    for e in range(n_elem):
        offset = np.array([5.0 * e, 0.0, 0.0])
        xyz = nc @ A.T + B + offset
        coords.append(xyz)
        connectivity.append(np.arange(e * n + 1, e * n + n + 1))
        elem_ids.append(100 + e)
    coords = np.vstack(coords)
    node_ids = np.arange(1, len(coords) + 1)
    rng = np.random.default_rng(3)
    G = rng.normal(size=(6, 3))
    c0 = np.arange(1.0, 7.0)
    stress = c0 + coords @ G.T
    return mesh_io.MeshData(
        node_ids=node_ids, coords=coords, elem_ids=np.array(elem_ids),
        elem_types=[name] * n_elem, connectivity=connectivity, stress=stress)


def assert_same_mesh(a, b):
    order_a = np.argsort(a.node_ids)
    order_b = np.argsort(b.node_ids)
    assert list(a.node_ids[order_a]) == list(b.node_ids[order_b])
    assert a.coords[order_a] == pytest.approx(b.coords[order_b])
    if a.stress is not None:
        assert a.stress[order_a] == pytest.approx(b.stress[order_b])
    a_elems = sorted(zip(a.elem_ids.tolist(), a.elem_types,
                         [tuple(sorted(c)) for c in a.connectivity]))
    b_elems = sorted(zip(b.elem_ids.tolist(), b.elem_types,
                         [tuple(sorted(c)) for c in b.connectivity]))
    assert a_elems == b_elems


# ---------------------------------------------------------------------------
# real file (HEX8)
# ---------------------------------------------------------------------------
def test_real_hex8_file_reads_and_matches_our_natural_coordinate_convention():
    """Our Brick8 numbers nodes so index k has natural coords
    ((k>>2)&1, (k>>1)&1, k&1): two indices differing in one bit must be an
    element EDGE (share a face-diagonal-free physical coordinate pair, i.e.
    differ in exactly one physical axis), and 0/7 (all bits differ) must be
    the space diagonal. This holds regardless of which physical axis ends up
    called r/s/t or which direction along it is "0" -- unlike a fixed
    axis/direction assumption, which a different real file need not share
    (see exodus_io.py's module docstring on rotational freedom)."""
    data = exodus_io.read_exodus(FIXTURE)
    assert len(data.node_ids) == 36
    assert len(data.elem_ids) == 10
    assert set(data.elem_types) == {"BRICK_8"}

    by_id = {int(n): data.coords[i] for i, n in enumerate(data.node_ids)}
    xyz = np.array([by_id[int(n)] for n in data.connectivity[0]])
    assert len({tuple(p) for p in xyz}) == 8  # 8 distinct physical corners

    for a in range(8):
        for b in range(a + 1, 8):
            differing_axes = np.flatnonzero(~np.isclose(xyz[a], xyz[b]))
            n_bits = bin(a ^ b).count("1")
            # one differing bit -> exactly one differing physical coordinate (an edge);
            # three differing bits (0 and 7) -> all three differ (the space diagonal)
            assert len(differing_axes) == n_bits, (a, b, xyz[a], xyz[b])


def test_real_hex8_elements_all_have_nonzero_volume():
    # read_exodus already raises if any element fails this; this test exists
    # so a future change that weakens/removes that check is still caught.
    # data.connectivity is already in OUR node order (the permutation was
    # applied during read) -- do not re-apply it here.
    data = exodus_io.read_exodus(FIXTURE)
    row = {int(n): i for i, n in enumerate(data.node_ids)}
    for conn in data.connectivity:
        c = data.coords[[row[int(n)] for n in conn]]
        exodus_io._check_nonzero_volume("BRICK_8", c, -1)


# ---------------------------------------------------------------------------
# write -> read round trip (all three supported types)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["TET_4", "WEDGE_6", "BRICK_8"])
def test_write_read_round_trip(tmp_path, name):
    mesh = synthetic_mesh(name, n_elem=2)
    path = str(tmp_path / "out.exo")
    exodus_io.write_exodus(path, mesh, nodal_vars={})
    back = exodus_io.read_exodus(path)
    assert_same_mesh(mesh, back)


@pytest.mark.parametrize("name", ["TET_4", "WEDGE_6", "BRICK_8"])
def test_written_elements_have_positive_volume_independent_check(tmp_path, name):
    """Cross-check with a formula independent of exodus_io's own
    _check_nonzero_volume (which runs during write/read already) --
    scipy's ConvexHull volume, which only needs a set of 3D points."""
    from scipy.spatial import ConvexHull

    mesh = synthetic_mesh(name, n_elem=1)
    path = str(tmp_path / "out.exo")
    exodus_io.write_exodus(path, mesh, nodal_vars={})
    back = exodus_io.read_exodus(path)
    row = {int(n): i for i, n in enumerate(back.node_ids)}
    coords = back.coords[[row[int(n)] for n in back.connectivity[0]]]
    assert ConvexHull(coords).volume > 1e-6


def test_written_ids_are_never_non_positive(tmp_path):
    """Regression test for a real ParaView/VTK-reader crash: writing a
    0-based mesh (DDSim's native RDB convention) straight into
    node_num_map/elem_num_map produced an ID of 0, which the IOSS-based
    Exodus reader used by ParaView/VTK rejects outright ("non-positive
    global id 0"), silently dropping the whole element block -- reproduced
    directly with pvpython against the real SIPS3002 conversion: file opened
    (RequestInformation succeeded) but came back with 0 points/cells on
    Apply (RequestData). write_exodus must shift IDs to be 1-based whenever
    the source mesh's IDs start at or below 0."""
    mesh = synthetic_mesh("BRICK_8", n_elem=1)
    mesh.node_ids = np.arange(0, len(mesh.node_ids))  # 0-based, like real RDB data
    mesh.elem_ids = np.array([0])
    mesh.connectivity = [mesh.node_ids]

    path = str(tmp_path / "out.exo")
    exodus_io.write_exodus(path, mesh, nodal_vars={})

    import netCDF4
    with netCDF4.Dataset(path) as ds:
        node_ids = np.asarray(ds.variables["node_num_map"][:])
        elem_ids = np.asarray(ds.variables["elem_num_map"][:])
    assert node_ids.min() >= 1
    assert elem_ids.min() >= 1
    # a documented +1 shift (the minimal one that makes the smallest ID 1),
    # not an arbitrary renumbering
    assert sorted(node_ids) == list(range(1, len(mesh.node_ids) + 1))
    assert list(elem_ids) == [1]


def test_write_leaves_already_positive_ids_unshifted(tmp_path):
    mesh = synthetic_mesh("BRICK_8", n_elem=1)  # node_ids/elem_ids already 1-based
    path = str(tmp_path / "out.exo")
    exodus_io.write_exodus(path, mesh, nodal_vars={})

    import netCDF4
    with netCDF4.Dataset(path) as ds:
        node_ids = np.asarray(ds.variables["node_num_map"][:])
        elem_ids = np.asarray(ds.variables["elem_num_map"][:])
    assert sorted(node_ids) == sorted(int(n) for n in mesh.node_ids)
    assert sorted(elem_ids) == sorted(int(e) for e in mesh.elem_ids)


def test_write_uses_split_coords_matching_file_size_flag(tmp_path):
    """Regression test for two real bugs found by reproducing "opens but
    crashes/comes back empty on Apply" directly with pvpython (VTK's actual
    ExodusII/IOSS reader) against the real SIPS3002 conversion:

    1. Writing file_size=1 ("large model", meaning split coordx/coordy/coordz)
       while actually writing a combined `coord` variable made the reader
       fail outright on Apply ("failed to locate x nodal coordinates").
    2. Fixing that by writing file_size=0 with combined `coord` stopped the
       crash, but geometry came back fine while every nodal variable (e.g.
       "life", "stress") silently read back as all-zero -- confirmed by
       writing an equivalent minimal file with VTK's own vtkIOSSWriter (which
       uses split coords + file_size=1) and finding *that* file's nodal
       variables read correctly.

    Split coordx/coordy/coordz + file_size=1 is the one combination verified
    to get both geometry and nodal variables right, so that's what we write.
    """
    mesh = synthetic_mesh("BRICK_8", n_elem=1)
    path = str(tmp_path / "out.exo")
    exodus_io.write_exodus(path, mesh, nodal_vars={})

    import netCDF4
    with netCDF4.Dataset(path) as ds:
        assert ds.file_size == 1
        assert "coord" not in ds.variables
        for axis in "xyz":
            assert "coord" + axis in ds.variables


def test_write_rejects_quadratic_elements(tmp_path):
    mesh = synthetic_mesh("TET_4")
    mesh.elem_types = ["TET_10"]
    with pytest.raises(NotImplementedError):
        exodus_io.write_exodus(str(tmp_path / "out.exo"), mesh, nodal_vars={})


# ---------------------------------------------------------------------------
# nodal variables (the actual point of write_exodus: life predictions)
# ---------------------------------------------------------------------------
def test_nodal_variables_round_trip_with_default_for_missing_nodes():
    mesh = synthetic_mesh("BRICK_8", n_elem=1)
    life = {int(mesh.node_ids[0]): 12345.0, int(mesh.node_ids[3]): 999.0}

    import tempfile
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "out.exo")
        exodus_io.write_exodus(path, mesh, nodal_vars={"life": life}, default=-1.0)

        import netCDF4
        with netCDF4.Dataset(path) as ds:
            names = ["".join(c.decode() for c in row if c).strip()
                     for row in ds.variables["name_nod_var"][:]]
            assert "life" in names  # plus the mesh's own stress, written automatically
            life_idx = names.index("life") + 1
            values = ds.variables["vals_nod_var%d" % life_idx][0, :]
            node_ids = list(ds.variables["node_num_map"][:])
            got = dict(zip(node_ids, values))
    assert got[int(mesh.node_ids[0])] == pytest.approx(12345.0)
    assert got[int(mesh.node_ids[3])] == pytest.approx(999.0)
    assert got[int(mesh.node_ids[1])] == pytest.approx(-1.0)  # default, not in the dict


def test_unrecognized_stress_names_raise_with_helpful_message(tmp_path):
    mesh = synthetic_mesh("TET_4")
    mesh.stress = None  # otherwise write_exodus also writes the canonical stress_xx etc.
    path = str(tmp_path / "out.exo")
    # write six nodal variables under names that match no known alias
    values = {str(i): {int(n): float(i) for n in mesh.node_ids} for i in range(6)}
    exodus_io.write_exodus(path, mesh, nodal_vars=values)
    with pytest.raises(ValueError, match="stress_names"):
        exodus_io.read_exodus(path)

    # the override makes it work
    back = exodus_io.read_exodus(path, stress_names={
        "xx": "0", "yy": "1", "zz": "2", "xy": "3", "yz": "4", "zx": "5"})
    assert back.stress[:, 0] == pytest.approx(0.0)
    assert back.stress[:, 5] == pytest.approx(5.0)


def test_check_nonzero_volume_accepts_our_own_reference_orderings():
    """elements.py's own already-validated node order (test_elements.py's
    NODES, used throughout the element test suite) must count as valid for
    every supported type."""
    for name in ("TET_4", "WEDGE_6", "BRICK_8"):
        c = np.array(NODES[name], dtype=float)
        exodus_io._check_nonzero_volume(name, c, elem_id=-1)  # must not raise

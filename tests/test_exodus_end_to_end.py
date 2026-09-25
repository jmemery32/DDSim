"""End-to-end check: a real DDSim crack-growth run, on a mesh built entirely
of Exodus-supported linear elements, written out and read back as Exodus.

example1 and every SIF_Verification fixture use quadratic elements (TET_10 /
BRICK_20), which are outside this pass's scope (see exodus_io.py's module
docstring) -- so this uses a small synthetic single-BRICK_8-element cube with
the same physical setup as example1 (20x20x20, uniform syy = 3.4, the same
material) instead, reusing example1.par's material directly.
"""
import os

import numpy as np
import pytest

from ddsim import DamMo, MeshTools, Parameters, Vec3D, exodus_io
from ddsim import mesh_io

PAR_DIR = os.path.join(os.path.dirname(__file__), "..", "examples", "example1") + os.sep


def cube_model():
    """One BRICK_8 element spanning [-10, 10]^3, uniform stress syy = 3.4 --
    physically identical to example1, just one linear element instead of 12
    quadratic tets."""
    corners = np.array([[x, y, z] for x in (-10, 10) for y in (-10, 10) for z in (-10, 10)])
    # elements.Brick8 node order: index k has natural coords ((k>>2)&1,(k>>1)&1,k&1)
    order = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
    coords = np.array([[-10 + 20 * r, -10 + 20 * s, -10 + 20 * t] for r, s, t in order], dtype=float)
    node_ids = np.arange(len(coords))
    stress = np.zeros((len(coords), 6))
    stress[:, 1] = 3.4  # syy, matching example1.par's setup exactly
    data = mesh_io.MeshData(
        node_ids=node_ids, coords=coords, elem_ids=np.array([0]),
        elem_types=["BRICK_8"], connectivity=[node_ids], stress=stress)
    return MeshTools.MeshTools(data=data)


def test_life_prediction_round_trips_through_exodus(tmp_path):
    model = cube_model()
    model.SetPointInsideTolerance(1.0e-7)
    params = Parameters.Parameters("example1", PAR_DIR)
    cracks = DamMo.DamModel(model, model.GetNodeList(), 0, 0, "0", 0, params, None, None, None)

    # a corner node: quarter-ellipse crack, same physical setup validated
    # (for the quadratic-element example1) in test_end_to_end.py
    doid = 0
    xyz = model.GetNodeInfo(doid)[0]
    cracks.AddDamOro(doid)
    cracks.AddFDam(doid, xyz, 0.01, 0.01, params.material, 0, 0)
    life = cracks.SimDamGrowth(doid, params, "err", 100.0, params.material[16])
    assert life > 0

    path = str(tmp_path / "out.exo")
    cracks.ToExodusFile(path)

    back = exodus_io.read_exodus(path)
    assert set(back.elem_types) == {"BRICK_8"}

    import netCDF4
    with netCDF4.Dataset(path) as ds:
        names = ["".join(c.decode() for c in row if c).strip()
                 for row in ds.variables["name_nod_var"][:]]
        life_idx = names.index("life") + 1
        vals = ds.variables["vals_nod_var%d" % life_idx][0, :]
        node_ids = list(ds.variables["node_num_map"][:])
    got = dict(zip(node_ids, vals))

    assert got[doid] == pytest.approx(life)
    # every other node was never run at all -> NaN (ParaView shows it as
    # masked/blank), distinct from cracks.HighLife (ran, never failed)
    others = [v for n, v in got.items() if n != doid]
    assert all(np.isnan(v) for v in others)

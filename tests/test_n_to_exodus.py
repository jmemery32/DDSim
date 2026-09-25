"""n_to_exodus: aggregate a real .N life-prediction result file onto its RDB
mesh as a nodal variable, without running any crack growth. Uses the same
hand-written tiny linear-element RDB fixture as test_rdb_to_exodus.py (a
one-BRICK_8 cube) since example1/SIF_Verification use quadratic elements.
"""
import numpy as np
import pytest

from ddsim import exodus_io
from ddsim.tools import n_to_exodus
from test_rdb_to_exodus import write_one_brick_rdb


def write_n_file(path, rows):
    with open(path, "w") as f:
        for doid, rid, n in rows:
            f.write("%d %d %d \n" % (doid, rid, n))


def test_read_n_file_groups_by_node(tmp_path):
    path = str(tmp_path / "results.N")
    write_n_file(path, [(0, 1, 100), (0, 2, 300), (5, -1, 100000)])
    per_node = n_to_exodus.read_n_file(path)
    assert per_node == {0: [100.0, 300.0], 5: [100000.0]}


def test_node_means_matches_plain_arithmetic_mean():
    per_node = {0: [10.0, 20.0, 30.0], 1: [100000.0]}
    means = n_to_exodus.node_means(per_node)
    assert means[0] == pytest.approx(20.0)
    assert means[1] == pytest.approx(100000.0)  # single-row "sentinel" node


def test_convert_writes_per_node_mean_life_and_nan_for_missing_nodes(tmp_path):
    rdb_base = str(tmp_path / "cube")
    write_one_brick_rdb(rdb_base)

    n_path = str(tmp_path / "results.N")
    # node 0: three particles -> mean 20000; node 3: one deterministic/no-particle
    # sentinel row (rid -1); nodes 1,2,4-7 never appear -> NaN.
    write_n_file(n_path, [
        (0, 101, 10000), (0, 202, 20000), (0, 303, 30000),
        (3, -1, 100000),
    ])

    out = str(tmp_path / "life.exo")
    mesh, values = n_to_exodus.convert(rdb_base, n_path, out, var_name="life")
    assert values[0] == pytest.approx(20000.0)
    assert values[3] == pytest.approx(100000.0)

    back = exodus_io.read_exodus(out)
    assert len(back.node_ids) == 8

    import netCDF4
    with netCDF4.Dataset(out) as ds:
        names = ["".join(c.decode() for c in row if c).strip()
                 for row in ds.variables["name_nod_var"][:]]
        life_idx = names.index("life") + 1
        life_vals = ds.variables["vals_nod_var%d" % life_idx][0, :]
        node_ids = list(ds.variables["node_num_map"][:])
    got = dict(zip(node_ids, life_vals))

    assert got[1] == pytest.approx(20000.0)  # doid 0 -> node_num_map shifts +1
    assert got[4] == pytest.approx(100000.0)  # doid 3 -> 4
    for doid in (1, 2, 4, 5, 6, 7):
        assert np.isnan(got[doid + 1])


def test_main_requires_three_or_four_arguments(monkeypatch, capsys):
    import sys
    monkeypatch.setattr(sys, "argv", ["n_to_exodus", "only_one_arg"])
    try:
        n_to_exodus.main()
        assert False, "should have exited"
    except SystemExit as e:
        assert e.code == 1
    assert "usage" in capsys.readouterr().out

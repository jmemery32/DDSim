"""exodus_to_n: the reverse of n_to_exodus -- pull a nodal variable back out
of an Exodus file as a .N-style text file. Uses the same hand-written tiny
linear-element RDB fixture as test_rdb_to_exodus.py/test_n_to_exodus.py.
"""
import pytest

from ddsim import exodus_io
from ddsim.tools import exodus_to_n, n_to_exodus
from test_rdb_to_exodus import write_one_brick_rdb


def test_read_nodal_variable_matches_what_was_written(tmp_path):
    # write_one_brick_rdb's node ids are 0-based (0-7), so write_exodus
    # shifts every id by +1 on disk (Exodus's external IDs must be
    # strictly positive) -- read_nodal_variable undoes that shift (see its
    # own docstring and exodus_io.write_exodus's ddsim_node_offset
    # attribute), so the original doids come back unchanged here.
    rdb_base = str(tmp_path / "cube")
    write_one_brick_rdb(rdb_base)
    from ddsim import mesh_io
    mesh = mesh_io.read_rdb(rdb_base)

    out = str(tmp_path / "life.exo")
    exodus_io.write_exodus(out, mesh, {"life": {1: 10.0, 5: 20.0}},
                           default=float("nan"))

    got = exodus_io.read_nodal_variable(out, "life")
    assert got[1] == pytest.approx(10.0)
    assert got[5] == pytest.approx(20.0)
    assert set(got) == set(mesh.node_ids)
    assert got[0] != got[0]  # NaN, untouched


def test_read_nodal_variable_recovers_original_ids_even_when_shifted_on_disk(tmp_path):
    """Regression test for a real bug found building this tool: the real
    SIPS3002 mesh has a node ID of 0 (some other, unrelated node), so
    write_exodus's "make IDs positive" shift silently applies to every
    real doid too, not just an originally-0-based toy mesh -- confirmed by
    comparing a direct .N-file mean against the value this tool reported
    for the same doid, which were off by exactly one doid's worth. Checks
    the raw on-disk IDs ARE shifted (so this isn't just testing a no-op),
    then that the high-level reader still recovers the original ones."""
    rdb_base = str(tmp_path / "cube")
    write_one_brick_rdb(rdb_base)  # node ids 0-7 -- triggers the +1 shift
    from ddsim import mesh_io
    mesh = mesh_io.read_rdb(rdb_base)

    out = str(tmp_path / "life.exo")
    exodus_io.write_exodus(out, mesh, {"life": {3: 42.0}}, default=float("nan"))

    import netCDF4
    with netCDF4.Dataset(out) as ds:
        raw_ids = sorted(int(n) for n in ds.variables["node_num_map"][:])
        assert raw_ids == list(range(1, 9))  # shifted: 0-7 -> 1-8
        assert int(ds.ddsim_node_offset) == 1

    got = exodus_io.read_nodal_variable(out, "life")
    assert got[3] == pytest.approx(42.0)  # original doid, not shifted to 4


def test_read_nodal_variable_raises_on_unknown_name(tmp_path):
    rdb_base = str(tmp_path / "cube")
    write_one_brick_rdb(rdb_base)
    from ddsim import mesh_io
    mesh = mesh_io.read_rdb(rdb_base)
    out = str(tmp_path / "life.exo")
    exodus_io.write_exodus(out, mesh, {"life": {1: 10.0}})

    with pytest.raises(KeyError, match="no nodal variable"):
        exodus_io.read_nodal_variable(out, "not_a_real_variable")


def test_round_trips_through_n_to_exodus_and_back(tmp_path):
    """.N -> Exodus (n_to_exodus) -> .N (exodus_to_n): the mean life for
    each doid survives the round trip, original doid numbers intact
    (NaN/"never run" nodes dropped, same convention as a doid that was
    never in doid_list -- see module docstring)."""
    rdb_base = str(tmp_path / "cube")
    write_one_brick_rdb(rdb_base)

    n_path = str(tmp_path / "results.N")
    with open(n_path, "w") as f:
        f.write("1 101 10000 \n1 202 30000 \n5 -1 100000 \n")

    exo_path = str(tmp_path / "life.exo")
    n_to_exodus.convert(rdb_base, n_path, exo_path, var_name="life")

    roundtrip_path = str(tmp_path / "roundtrip.N")
    values = exodus_to_n.convert(exo_path, roundtrip_path, var_name="life")
    assert values[1] == pytest.approx(20000.0)  # mean of 10000, 30000
    assert values[5] == pytest.approx(100000.0)

    lines = open(roundtrip_path).read().splitlines()
    parsed = {}
    for line in lines:
        doid, rid, n = line.split()
        assert int(rid) == -1  # the only per-doid granularity this format has
        parsed[int(doid)] = float(n)
    assert parsed == {1: pytest.approx(20000.0), 5: pytest.approx(100000.0)}
    # nodes that never appeared in the original .N (NaN) are dropped, not
    # written as a literal row -- there are 8 nodes total, only 2 appear.
    assert len(lines) == 2


def test_main_requires_two_or_three_arguments(monkeypatch, capsys):
    import sys
    monkeypatch.setattr(sys, "argv", ["exodus_to_n", "only_one_arg"])
    try:
        exodus_to_n.main()
        assert False, "should have exited"
    except SystemExit as e:
        assert e.code == 1
    out = capsys.readouterr().out
    assert "sage" in out  # "Usage"/"usage", from the module docstring

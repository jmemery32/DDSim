"""rdb_to_exodus: a real RDB mesh on disk (not a synthetic MeshData) converted
to Exodus. Writes a tiny linear-element RDB fixture by hand -- example1 and
the SIF_Verification fixtures all use quadratic elements, which write_exodus
doesn't support (see exodus_io.py), so they can't exercise this path.
"""
import os

from ddsim.tools import rdb_to_exodus


def write_one_brick_rdb(base):
    # elements.Brick8 node order: index k has natural coords ((k>>2)&1,(k>>1)&1,k&1)
    order = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
    with open(base + ".nod", "w") as f:
        for i, (r, s, t) in enumerate(order):
            f.write("%d %g %g %g\n" % (i, -10 + 20 * r, -10 + 20 * s, -10 + 20 * t))
    with open(base + ".con", "w") as f:
        f.write("0 1 0 8 " + " ".join(str(i) for i in range(8)) + "\n")
    with open(base + ".smp", "w") as f:
        f.write("1 Brick1IShapeFunc\n")
    with open(base + ".sig", "w") as f:
        for i in range(8):
            f.write("%d 0.0 3.4 0.0 0.0 0.0 0.0\n" % i)  # uniform syy = 3.4, like example1


def test_convert_writes_a_readable_exodus_file(tmp_path):
    base = str(tmp_path / "cube")
    write_one_brick_rdb(base)
    out = str(tmp_path / "cube.exo")

    data = rdb_to_exodus.convert(base, out)
    assert len(data.node_ids) == 8
    assert os.path.exists(out)

    from ddsim import exodus_io
    back = exodus_io.read_exodus(out)
    assert len(back.node_ids) == 8
    assert set(back.elem_types) == {"BRICK_8"}
    assert (back.stress[:, 1] == 3.4).all()  # syy round-tripped


def test_main_requires_two_arguments(monkeypatch, capsys):
    import sys
    monkeypatch.setattr(sys, "argv", ["rdb_to_exodus", "only_one_arg"])
    try:
        rdb_to_exodus.main()
        assert False, "should have exited"
    except SystemExit as e:
        assert e.code == 1
    assert "usage" in capsys.readouterr().out

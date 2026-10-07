"""Tests for -crack_path: writing a doid's full crack-growth history as a
legacy VTK PolyData file, for visualizing the predicted crack path directly
in ParaView (see DamMo.DamModel.WriteCrackPathVTK, docs/PORTING_NOTES.md).
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

import pytest

from ddsim import DamMo, Parameters, VarAmplitude
from test_exodus_end_to_end import cube_model, PAR_DIR

ROOT = os.path.join(os.path.dirname(__file__), "..")
EXAMPLE = os.path.join(ROOT, "examples", "example1")


def parse_legacy_vtk_polydata(path):
    """Minimal hand-rolled parser for exactly what WriteCrackPathVTK writes:
    POINTS, LINES, and three CELL_DATA SCALARS arrays (N, crack_type, step).
    Deliberately independent of any VTK/ParaView library."""
    lines = open(path).read().splitlines()
    i = 0
    n_points = n_lines = None
    points = []
    cells = []
    cell_data = {}

    while i < len(lines):
        parts = lines[i].split()
        if parts[:1] == ["POINTS"]:
            n_points = int(parts[1])
            for j in range(n_points):
                points.append(tuple(float(x) for x in lines[i + 1 + j].split()))
            i += 1 + n_points
        elif parts[:1] == ["LINES"]:
            n_lines = int(parts[1])
            for j in range(n_lines):
                vals = [int(x) for x in lines[i + 1 + j].split()]
                cells.append(vals[1:1 + vals[0]])
            i += 1 + n_lines
        elif parts[:1] == ["SCALARS"]:
            name = parts[1]
            n_cells = len(cells)
            values = [float(x) for x in lines[i + 2:i + 2 + n_cells]]
            cell_data[name] = values
            i += 2 + n_cells
        else:
            i += 1
    assert n_points is not None and n_lines is not None
    return points, cells, cell_data


def run_with_crack_path(tmp_path, nodes, scale, extra_args=()):
    work = tmp_path / "example1"
    shutil.copytree(EXAMPLE, work)
    out = work / "crack_path.vtk"
    proc = subprocess.run(
        [sys.executable, "-m", "ddsim.DDSim", "-base", "example1",
         "-conpath", "./", "-parpath", "./", "-v", "-doid_list", nodes,
         "-scale", scale, "-crack_path", str(out), *extra_args],
        cwd=work, capture_output=True, text=True, timeout=600)
    return proc, out


def test_crack_path_vtk_matches_reported_life(tmp_path):
    proc, out = run_with_crack_path(tmp_path, "10", "100")
    assert proc.returncode == 0, proc.stderr[-2000:]
    m = re.search(r"Life is:\s+(\S+)", proc.stdout)
    life = float(m.group(1))

    points, cells, cell_data = parse_legacy_vtk_polydata(out)
    assert len(cells) > 1  # more than just the initial state
    assert cell_data["N"][-1] == pytest.approx(life, rel=1e-8)
    assert cell_data["N"] == sorted(cell_data["N"])  # monotonically increasing
    assert cell_data["step"] == list(range(len(cells)))
    assert set(cell_data["crack_type"]) <= {0.0, 1.0, 2.0}

    # every line indexes real points, and each index list is internally
    # distinct except the deliberate Fellipse closing repeat (first==last)
    for cell in cells:
        assert all(0 <= idx < len(points) for idx in cell)
        interior = cell[:-1] if cell[0] == cell[-1] else cell
        assert len(set(interior)) == len(interior)


def test_crack_path_fellipse_steps_close_the_loop(tmp_path):
    """Fellipse (fully embedded) is the only crack type whose polyline
    should repeat its first point index to close the loop visually --
    Hellipse/Qellipse are genuinely open arcs ending at the free surface."""
    _, out = run_with_crack_path(tmp_path, "10", "100")
    _, cells, cell_data = parse_legacy_vtk_polydata(out)
    for cell, ctype in zip(cells, cell_data["crack_type"]):
        if ctype == 0.0:  # Fellipse
            assert cell[0] == cell[-1]
        else:  # Hellipse/Qellipse
            assert cell[0] != cell[-1]


def test_crack_path_errors_clearly_under_dash_j(tmp_path):
    """-j workers strip DamEl down to just the final state (see
    _StrippedDamEl) -- the full history this needs is gone by the time
    results reach the parent, so this must fail loudly, not silently write
    a useless file."""
    proc, out = run_with_crack_path(tmp_path, "10,9", "100", extra_args=["-j", "2"])
    assert proc.returncode != 0
    assert "crack-growth history isn't available" in proc.stderr
    assert not out.exists()


# ---------------------------------------------------------------------------
# state_N: the absolute-cumulative-cycle-count bookkeeping -crack_path reads
# (see docs/PORTING_NOTES.md, "Bug found later: -crack_path's cumulative N
# was wrong for VarAmp"). Exercises VarAmp directly (not via the CLI) since
# that's the path the real bug was in; GrowDam/RK5 is already covered above
# via the full subprocess tests.
# ---------------------------------------------------------------------------
def test_state_n_is_index_aligned_and_matches_reported_life():
    model = cube_model()
    model.SetPointInsideTolerance(1.0e-7)
    params = Parameters.Parameters("example1", PAR_DIR)
    cracks = DamMo.DamModel(model, model.GetNodeList(), 0, False, "0",
                            False, params, False, False, "RK5")
    doid = 0
    xyz = model.GetNodeInfo(doid)[0]
    cracks.AddDamOro(doid)
    cracks.AddFDam(doid, xyz, 0.001, 0.001, params.material, 0, False)

    spec_text = "\n".join("%f %f" % (0.05 * (i % 4), 1.0) for i in range(20)) + "\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".val", delete=False) as tf:
        tf.write(spec_text)
        val_path = tf.name
    try:
        spec = VarAmplitude.Spectrum(val_path)
        N, _fin = cracks.VarAmp(doid, params.Max_crack_size, 50000, spec,
                                0.0, 100.0, params.r)
    finally:
        os.unlink(val_path)

    dam_els = cracks.DamOro[doid].DamEl
    assert len(dam_els) >= 2, "expected a Fellipse->Qellipse transition here"

    prev_last_state_n = None
    for stage in dam_els:
        # the core invariant WriteCrackPathVTK depends on: one state_N
        # entry per recorded state, never one per self.dN entry (which,
        # for VarAmp, can be longer -- see the stall-cycle discussion in
        # PORTING_NOTES.md).
        assert len(stage.state_N) == len(stage.a[0][0])
        assert stage.state_N == sorted(stage.state_N)  # non-decreasing
        if prev_last_state_n is not None:
            # regression check for the transition-carryover bug: a new
            # stage's first state_N must inherit the previous stage's last
            # one (the same physical state), not reset to 0.0.
            assert stage.state_N[0] == pytest.approx(prev_last_state_n)
        prev_last_state_n = stage.state_N[-1]

    assert dam_els[-1].state_N[-1] == pytest.approx(N)

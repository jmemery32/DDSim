"""Golden end-to-end regression: run the real ddsim.DDSim driver on example1.

The numbers below are the life predictions of the pure-Python (pre-numba)
implementation for nodes 10 and 0, which the compiled kernels reproduce to the
last digit.  Nodes 12 and 9 were recorded from the compiled code.  They guard
every later optimisation: performance work must not change these.

Scales are deliberately high (K scale factor on the uniform 3.4 stress) so the
cracks reach unstable growth quickly:  node 0 = corner (quarter-ellipse crack),
node 10 = face centre (half-ellipse), nodes 12 / 9 = edge nodes.

Re-recorded 2026 (see docs/PORTING_NOTES.md, "RK5 near-instability overshoot")
after fixing a real bug where GrowDam's RK5 path could take a single step
large enough to blow past the stable-growth region entirely. The fix is a
direct cap on relative crack growth per step (Parameters.max_growth_fraction,
default 10%); every value below reflects the fixed, slower-stepping behavior
near instability, confirmed sane via independent cross-checks against plain
forward-Euler and a cycle-by-cycle VarAmp integration on the same geometry.
"""
import os
import re
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
EXAMPLE = os.path.join(ROOT, "examples", "example1")

# (nodes, scale) -> {doid: (life, will_grow)}
GOLDEN = {
    ("10", "100"): {10: (27.02528469667873, 2)},
    ("0,12,9", "60"): {
        0: (406.6135133294466, 4),
        12: (117.7760415435534, 4),
        9: (107.23339738673845, 2),
    },
}
LINE = re.compile(r"doid:\s+(\d+)\s+Life is:\s+(\S+)\s+WillGrow:\s+(\d+)")


def run_driver(tmp_path, nodes, scale, extra_args=()):
    work = tmp_path / "example1"
    shutil.copytree(EXAMPLE, work)
    proc = subprocess.run(
        [sys.executable, "-m", "ddsim.DDSim", "-base", "example1",
         "-conpath", "./", "-parpath", "./", "-v", "-doid_list", nodes, "-scale", scale,
         *extra_args],
        cwd=work, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stderr[-2000:]
    return {int(d): (float(life), int(wg)) for d, life, wg in LINE.findall(proc.stdout)}


@pytest.mark.parametrize("nodes,scale", list(GOLDEN))
def test_life_predictions_match_golden(tmp_path, nodes, scale):
    got = run_driver(tmp_path, nodes, scale)
    expected = GOLDEN[(nodes, scale)]
    assert set(got) == set(expected)
    for doid, (life, will_grow) in expected.items():
        assert got[doid][1] == will_grow
        assert got[doid][0] == pytest.approx(life, rel=1e-9)


def test_cube_symmetry_all_eight_corners_have_the_same_life(tmp_path):
    """The cube and its uniform stress are symmetric, so the eight corner nodes must
    give the same life whatever their orientation relative to the load.  This
    catches direction-dependent bugs (surface normals, crack frames, mesh
    orientation) that a single-node golden value cannot."""
    got = run_driver(tmp_path, "0,1,2,3,4,5,6,7", "100")
    assert sorted(got) == list(range(8))
    lives = [got[d][0] for d in range(8)]
    assert all(got[d][1] == 4 for d in range(8))
    assert lives == pytest.approx([lives[0]] * 8, rel=1e-8)
    assert lives[0] == pytest.approx(79.14869788864108, rel=1e-8)


def test_interior_node_with_uniaxial_stress_does_not_recurse_forever(tmp_path):
    """Node 8 (cube centre); under uniaxial stress sigma_2 = sigma_3, so
    re-orienting the crack never helps.  The 2007 code died with a RecursionError.

    Before the 2026 RK5 near-instability fix (see docs/PORTING_NOTES.md), this
    node's crack took one oversized RK5 step that overshot clean past the
    stable-growth region and outgrew the whole cube (WillGrow 4, "net
    fracture") before the code could ever classify it as unstable. The fixed,
    growth-rate-capped stepping now catches the same crack while it's still
    inside the cube, correctly classifying it as unstable growth (WillGrow 2)
    at a much smaller, physically sane size -- this is the better answer, not
    a regression; the RecursionError this test was written to guard against
    is unrelated and still guarded by the subprocess completing at all."""
    got = run_driver(tmp_path, "8", "100")
    assert got[8][1] == 2
    assert got[8][0] == pytest.approx(33.42699710068892, rel=1e-8)

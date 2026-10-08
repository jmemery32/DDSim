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

# These golden values were captured on one specific machine/library-version
# combination. Everything here sits right at (nodes 10, 8) or near (the cube
# corners, node 9) the RK5 near-instability region this port's GrowDam fix
# deliberately changed (see docs/PORTING_NOTES.md) -- genuinely chaotic:
# different numpy/scipy/numba builds (different BLAS, different numba JIT
# codegen) can nudge which side of a step-size decision a run lands on,
# producing a real but small (<1%, observed up to ~0.6%) drift in the final
# answer, with no code bug involved. CROSS_ENV_TOL absorbs that so CI running
# multiple OS/Python versions doesn't flag it, while still catching any real
# regression -- which this session's actual bugs moved these values by 10x-
# 1000x, not fractions of a percent.
CROSS_ENV_TOL = 1e-2

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
        assert got[doid][0] == pytest.approx(life, rel=CROSS_ENV_TOL)


def test_cube_symmetry_all_eight_corners_have_the_same_life(tmp_path):
    """The cube and its uniform stress are symmetric, so the eight corner nodes must
    give the same life whatever their orientation relative to the load.  This
    catches direction-dependent bugs (surface normals, crack frames, mesh
    orientation) that a single-node golden value cannot."""
    got = run_driver(tmp_path, "0,1,2,3,4,5,6,7", "100")
    assert sorted(got) == list(range(8))
    lives = [got[d][0] for d in range(8)]
    assert all(got[d][1] == 4 for d in range(8))
    assert lives == pytest.approx([lives[0]] * 8, rel=1e-8)  # live self-comparison: stays tight
    assert lives[0] == pytest.approx(79.14869788864108, rel=CROSS_ENV_TOL)


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
    assert got[8][0] == pytest.approx(33.42699710068892, rel=CROSS_ENV_TOL)


# ---------------------------------------------------------------------------
# -ai/-bi: command-line override of the deterministic initial flaw size
# ---------------------------------------------------------------------------
def test_ai_override_changes_initial_crack_size(tmp_path):
    """-ai overrides the .par file's a_b (example1.par's default is 0.01)."""
    default = run_driver(tmp_path / "default", "10", "100")
    overridden = run_driver(tmp_path / "overridden", "10", "100", extra_args=["-ai", "0.05"])
    assert default[10] != overridden[10]
    assert overridden[10] == (13.0, 4)


def test_bi_overrides_second_dimension_independently(tmp_path):
    """-bi sets the second crack dimension separately from -ai (an
    asymmetric initial flaw); omitting it (previous test) uses -ai for both."""
    got = run_driver(tmp_path, "10", "100", extra_args=["-ai", "0.05", "-bi", "0.08"])
    assert got[10] == (11.0, 2)


def test_ai_override_matches_across_j(tmp_path):
    """The -ai override must reach -j workers too -- they re-read the .par
    file independently (see ddsim.parallel), so the parent's own override of
    its in-memory Parameters.a_b doesn't reach them for free."""
    serial = run_driver(tmp_path / "serial", "10,9", "100", extra_args=["-ai", "0.05"])
    parallel = run_driver(tmp_path / "parallel", "10,9", "100", extra_args=["-ai", "0.05", "-j", "2"])
    assert serial == parallel == {10: (13.0, 4), 9: (12.0, 2)}


def test_ai_ignored_with_warning_for_monte_runs(tmp_path):
    """-ai only means something for monte=0 (a single literal flaw size);
    for monte=1/2, a_b sets the sampling distribution's shape instead, so the
    override is ignored with a printed warning rather than silently changing
    the distribution."""
    work = tmp_path / "example1"
    shutil.copytree(EXAMPLE, work)
    par = work / "example1.par"
    par.write_text(par.read_text().replace("monte\n0\n", "monte\n1\n"))
    # same doid set/scale/seed test_parallel.py's monte-carlo test already
    # confirms works for monte=1 -- a_b itself isn't exercised by this test,
    # just that -ai is ignored (with a warning) rather than changing anything.
    proc = subprocess.run(
        [sys.executable, "-m", "ddsim.DDSim", "-base", "example1",
         "-conpath", "./", "-parpath", "./", "-doid_list", "0,1,2,3,4,5,6,7",
         "-scale", "60", "-ai", "0.05", "-seed", "42"],
        cwd=work, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "WARNING" in proc.stdout and "-ai/-bi" in proc.stdout

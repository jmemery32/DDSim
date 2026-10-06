"""Tests for the -j multiprocessing rewrite (replacing the old Windows/MPI
cluster workflow -- see ddsim.parallel and docs/PORTING_NOTES.md) and the two
bug fixes that were prerequisites for it:

* DamMo.SimDamGrowth's .Life assignment (only affects deterministic runs).
* DDSim.Var_Amplitude's `verify` NameError and missing deterministic .Life.
"""
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from ddsim import DamMo, Parameters, parallel
from test_end_to_end import EXAMPLE, run_driver
from test_exodus_end_to_end import cube_model, PAR_DIR

ROOT = os.path.join(os.path.dirname(__file__), "..")


# ---------------------------------------------------------------------------
# partition_doids: in-memory equivalent of the old DDSim.DoParallel
# ---------------------------------------------------------------------------
def test_partition_doids_covers_every_doid_exactly_once():
    doids = list(range(37))
    chunks = parallel.partition_doids(doids, 5, seed=7)
    assert len(chunks) == 5
    flat = sorted(sum(chunks, []))
    assert flat == doids


def test_partition_doids_does_not_mutate_input():
    doids = list(range(10))
    original = list(doids)
    parallel.partition_doids(doids, 3, seed=1)
    assert doids == original


def test_partition_doids_matches_hand_computed_chunks():
    # Locks in the exact ported shuffle+chunk algorithm (same shuffle loop
    # and contiguous-slice split as the legacy DoParallel) for a fixed seed.
    chunks = parallel.partition_doids(list(range(10)), 3, seed=123)
    assert chunks == [sorted(c) for c in chunks]  # each chunk internally sorted
    assert [len(c) for c in chunks] == [3, 3, 4]  # 10 // 3 = 3, remainder to last chunk
    assert parallel.partition_doids(list(range(10)), 3, seed=123) == chunks  # deterministic


def test_partition_doids_handles_more_workers_than_doids():
    chunks = parallel.partition_doids([1, 2], 5, seed=1)
    assert sorted(sum(chunks, [])) == [1, 2]
    # some chunks may be empty; callers are expected to filter those out
    # (run_parallel clamps num_workers <= len(doid_list) before this point).


# ---------------------------------------------------------------------------
# StripDamElForTransport: the picklability fix
# ---------------------------------------------------------------------------
def test_strip_dam_el_for_transport_round_trips_through_pickle():
    import pickle

    model = cube_model()
    model.SetPointInsideTolerance(1.0e-7)
    params = Parameters.Parameters("example1", PAR_DIR)
    cracks = DamMo.DamModel(model, model.GetNodeList(), 0, 0, "0", 0, params, None, None, None)
    doid = 0
    xyz = model.GetNodeInfo(doid)[0]
    cracks.AddDamOro(doid)
    cracks.AddFDam(doid, xyz, 0.01, 0.01, params.material, 0, 0)
    life = cracks.SimDamGrowth(doid, params, "err", 100.0, params.material[16])

    rotation_before = cracks.DamOro[doid].DamEl[0].Rotation
    af_before, bf_before = cracks.DamOro[doid].DamEl[-1].GiveCurrent()

    # A live DamOro[doid] IS picklable now that _DamOroContainer is a
    # top-level class (that was the crash fix) -- but its DamEl holds a
    # direct reference to the whole shared mesh (via each DamEl's .model),
    # so it's pickled along too, ballooning the payload. That's what
    # StripDamElForTransport actually fixes: not an error, a size/cost
    # problem -- confirm the size drops by orders of magnitude, not just
    # that pickling succeeds. See docs/PORTING_NOTES.md.
    unstripped_size = len(pickle.dumps(cracks.DamOro[doid]))

    cracks.DamOro[doid].StripDamElForTransport()
    stripped_bytes = pickle.dumps(cracks.DamOro[doid])
    assert len(stripped_bytes) < unstripped_size / 10
    restored = pickle.loads(stripped_bytes)

    assert len(restored.DamEl) == 1
    assert restored.DamEl[0].Rotation == rotation_before
    assert restored.DamEl[-1].GiveCurrent() == (af_before, bf_before)
    assert restored.Life == pytest.approx(life)

    # what WriteRotations/WriteFinalAs actually read still matches
    import io
    cracks.DamOro[doid] = restored
    rot_file = io.StringIO()
    cracks.WriteRotations(rot_file, doid)
    af_file = io.StringIO()
    cracks.WriteFinalAs(af_file, doid)
    assert str(af_before) in af_file.getvalue()
    assert str(bf_before) in af_file.getvalue()


# ---------------------------------------------------------------------------
# SimDamGrowth .Life bug (deterministic-mode-only, original-2007 bug)
# ---------------------------------------------------------------------------
def test_simdamgrowth_sets_life_for_a_middle_of_the_pack_doid():
    """Regression test for the original-2007 bug: .Life used to only be set
    when a doid's computed N also happened to be a new model-wide running
    max/min. Process three doids (via three different initial crack sizes at
    the same physical point, since DamOro keys are just labels -- AddFDam
    takes xyz explicitly) through one DamModel, ordered so the third is
    neither a new high nor a new low."""
    model = cube_model()
    model.SetPointInsideTolerance(1.0e-7)
    params = Parameters.Parameters("example1", PAR_DIR)
    cracks = DamMo.DamModel(model, model.GetNodeList(), 0, 0, "0", 0, params, None, None, None)
    xyz = model.GetNodeInfo(0)[0]

    def run(doid, a):
        cracks.AddDamOro(doid)
        cracks.AddFDam(doid, xyz, a, a, params.material, 0, 0)
        return cracks.SimDamGrowth(doid, params, "err", 100.0, params.material[16])

    life_high = run(100, 0.005)  # smallest initial crack -> largest life -> new HighLife
    life_low = run(101, 0.05)    # largest initial crack -> smallest life -> new LowLife
    life_mid = run(102, 0.01)    # in between -> neither a new high nor a new low

    assert life_low < life_mid < life_high  # confirms the ordering this test relies on
    assert cracks.HighLife == pytest.approx(life_high)
    assert cracks.LowLife == pytest.approx(life_low)
    # the actual regression check: doid 102's own computed life, not -1
    assert cracks.DamOro[102].Life == pytest.approx(life_mid)
    assert cracks.DamOro[102].Life != -1


def test_refresh_life_bounds_recomputes_from_merged_damoro():
    model = cube_model()
    model.SetPointInsideTolerance(1.0e-7)
    params = Parameters.Parameters("example1", PAR_DIR)
    cracks = DamMo.DamModel(model, model.GetNodeList(), 0, 0, "0", 0, params, None, None, None)
    # simulate a post-merge DamModel: HighLife/LowLife still at their fresh
    # defaults (as the parent's would be after a -j run, since the parent
    # itself never calls SimDamGrowth), DamOro populated "as if" merged in.
    for doid, life in ((1, 50.0), (2, 10.0), (3, 30.0)):
        cracks.AddDamOro(doid)
        cracks.DamOro[doid].Life = life
    cracks.RefreshLifeBounds()
    assert cracks.HighLife == pytest.approx(50.0)
    assert cracks.LowLife == pytest.approx(10.0)


# ---------------------------------------------------------------------------
# End-to-end: -j reproduces the exact serial results
# ---------------------------------------------------------------------------
def test_j_matches_serial_golden_values(tmp_path):
    """The same doid set/scale that test_end_to_end.py's golden test uses,
    but with -j 3 (3 doids, 3 workers -- one doid per worker)."""
    got = run_driver(tmp_path, "0,12,9", "60", extra_args=["-j", "3"])
    expected = {
        0: (406.6135133294466, 4),
        12: (117.7760415435534, 4),
        9: (107.23339738673845, 2),
    }
    assert set(got) == set(expected)
    for doid, (life, will_grow) in expected.items():
        assert got[doid][1] == will_grow
        assert got[doid][0] == pytest.approx(life, rel=1e-9)


def test_j_matches_serial_with_more_doids_than_workers(tmp_path):
    """8 symmetric cube corners split across 3 workers -- some workers must
    handle multiple doids, exercising the chunking (not just 1 doid/worker)."""
    got = run_driver(tmp_path, "0,1,2,3,4,5,6,7", "100", extra_args=["-j", "3"])
    assert sorted(got) == list(range(8))
    assert all(got[d][1] == 4 for d in range(8))
    lives = [got[d][0] for d in range(8)]
    assert lives == pytest.approx([79.14869788864108] * 8, rel=1e-8)


def _run_saveall(tmp_path, name, nodes, scale, extra_args, patch_monte=False):
    work = tmp_path / name
    shutil.copytree(EXAMPLE, work)
    if patch_monte:
        par = work / "example1.par"
        par.write_text(par.read_text().replace("monte\n0\n", "monte\n1\n"))
    proc = subprocess.run(
        [sys.executable, "-m", "ddsim.DDSim", "-base", "example1",
         "-conpath", "./", "-parpath", "./", "-doid_list", nodes, "-scale", scale,
         "-sv", *extra_args],
        cwd=work, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stderr[-2000:]
    return work


def test_j_saveall_files_match_serial_deterministic(tmp_path):
    serial = _run_saveall(tmp_path, "serial", "0,1,2,3,4,5,6,7", "100", [])
    parallel_run = _run_saveall(tmp_path, "parallel", "0,1,2,3,4,5,6,7", "100", ["-j", "3"])
    for ext in (".N", ".ai", ".af", ".ori"):
        a = sorted((serial / ("example1" + ext)).read_text().splitlines())
        b = sorted((parallel_run / ("example1" + ext)).read_text().splitlines())
        assert a == b, ext


def test_j_saveall_files_match_serial_monte_carlo(tmp_path):
    """Same seed -> the random initial a's are generated once in the parent
    (before -j even branches), so results should be bit-identical to serial
    regardless of how doids are split across workers."""
    serial = _run_saveall(tmp_path, "serial", "0,1,2,3,4,5,6,7", "60",
                          ["-seed", "42"], patch_monte=True)
    parallel_run = _run_saveall(tmp_path, "parallel", "0,1,2,3,4,5,6,7", "60",
                                ["-seed", "42", "-j", "3"], patch_monte=True)
    a = sorted((serial / "example1.N").read_text().splitlines())
    b = sorted((parallel_run / "example1.N").read_text().splitlines())
    assert a == b


# ---------------------------------------------------------------------------
# Var_Amplitude bug fixes (verify NameError; missing deterministic .Life)
# ---------------------------------------------------------------------------
VAL_SPECTRUM = "\n".join("%f %f" % (0.1 * (i % 3), 1.0) for i in range(20)) + "\n"


def test_var_amplitude_deterministic_runs_without_crashing_and_sets_life(tmp_path):
    work = tmp_path / "example1"
    shutil.copytree(EXAMPLE, work)
    (work / "example1.val").write_text(VAL_SPECTRUM)
    proc = subprocess.run(
        [sys.executable, "-m", "ddsim.DDSim", "-base", "example1",
         "-conpath", "./", "-parpath", "./", "-doid_list", "0", "-VarAmp",
         "example1.val", "-sv"],
        cwd=work, capture_output=True, text=True, timeout=600)
    # regression check for the `verify` NameError: this used to crash
    # immediately on the first AddFDam call, every time -VarAmp was used.
    assert proc.returncode == 0, proc.stderr[-2000:]
    nfile = (work / "example1.N").read_text().split()
    # doid, rid, N triple; regression check for the missing .Life assignment
    # in the deterministic branch (used to always write -1 here).
    assert int(nfile[0]) == 0
    assert int(nfile[2]) != -1


def test_var_amplitude_monte_carlo_runs_without_crashing(tmp_path):
    """Covers the three Monte Carlo-branch AddFDam(...,verify) call sites the
    deterministic test above doesn't reach."""
    work = tmp_path / "example1"
    shutil.copytree(EXAMPLE, work)
    (work / "example1.val").write_text(VAL_SPECTRUM)
    par = work / "example1.par"
    par.write_text(par.read_text().replace("monte\n0\n", "monte\n1\n"))
    proc = subprocess.run(
        [sys.executable, "-m", "ddsim.DDSim", "-base", "example1",
         "-conpath", "./", "-parpath", "./", "-doid_list", "0", "-seed", "42",
         "-VarAmp", "example1.val", "-sv"],
        cwd=work, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stderr[-2000:]

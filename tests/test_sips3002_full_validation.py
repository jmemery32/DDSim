"""Validates the -j multiprocessing rewrite (and the SimDamGrowth/Var_Amplitude
bug fixes that were prerequisites for it) against the real captured 2007
particle-filter Monte Carlo results for the SIPS3002 open-hole model.

Two tiers, both using real, proprietary data kept OUTSIDE the repository
(see docs/PORTING_NOTES.md's "Multiprocessing" section for the full writeup
and the analysis that produced the thresholds hard-coded below):

* A fast (~20s) test that actually runs `ddsim -j 2` against the real
  particle-filter inputs for a small, representative doid sample (the known
  hot-spot cluster from Sec. 6.1's peak-stress cluster, plus several
  "sentinel"/no-particle-nearby nodes) and compares against the real
  captured `.N` file. Runs whenever `DDSIM_SIPS3002_DIR` is set (same env
  var `test_sips3002_sanity.py` uses).
* A slow, whole-model statistical check (all 63,974 nodes) against a
  *precomputed* `.N` file from a full `-S -sv -j <N>` run -- never runs the
  hour-long simulation itself. Opt in with `DDSIM_SIPS3002_FULL_N` pointing
  at that file's path; skipped otherwise.
"""
import os
import shutil
import subprocess
import sys

import pytest

from ddsim.tools import n_to_exodus

SIPS_DIR = os.environ.get("DDSIM_SIPS3002_DIR")
PARTICLES_DIR = os.path.join(SIPS_DIR, "10000_Particles") if SIPS_DIR else None
ORIGINAL_N = (os.path.join(PARTICLES_DIR, "3002_open_CA_10000Parts_070530.N")
              if PARTICLES_DIR else None)

pytestmark = pytest.mark.skipif(
    not SIPS_DIR or not ORIGINAL_N or not os.path.exists(ORIGINAL_N)
    or not os.path.exists(os.path.join(PARTICLES_DIR, "sips3002.par")),
    reason="set DDSIM_SIPS3002_DIR to a local copy of SIPS3002_open/ConstantAmplitude "
           "(needs 10000_Particles/{sips3002.par,.rnd,.map,"
           "3002_open_CA_10000Parts_070530.N} too) to run these")

# The known peak-stress hot-spot cluster (Sec. 6.1; see test_sips3002_sanity.py
# and PORTING_NOTES.md) -- the only nodes in the model with enough particles
# (thousands) to be affected by the input-data-provenance limitation
# documented in PORTING_NOTES.md (the .rnd/.map files on disk aren't a
# byte-exact match to whatever snapshot generated this historical result).
HOT_SPOT_DOIDS = [119032, 113831, 124125, 111308]

# Empirically established (see PORTING_NOTES.md) on the actual 172,601-node
# model, comparing a full run against 3002_open_CA_10000Parts_070530.N:
# 96.0% of all 63,974 nodes match to machine precision, 99.4% within 1%,
# ~100% within 10%, mean signed relative difference ~1.3e-6 (no bias).
# Thresholds below are deliberately looser than the observed numbers, so a
# minor future edit to the external data doesn't make this flaky, while
# still catching a real regression (e.g. a reintroduced order-dependence
# bug, which would blow well past all of these).
MIN_FRAC_MACHINE_PRECISION = 0.90
MIN_FRAC_WITHIN_1PCT = 0.98
MIN_FRAC_WITHIN_10PCT = 0.999
MAX_ABS_MEAN_SIGNED_RELDIFF = 0.01


def _read_original_n():
    """{doid: [(rid, N), ...]}, keeping rid (n_to_exodus.read_n_file drops
    it) since distinguishing a true rid=-1 sentinel from a real single
    particle matters here -- see PORTING_NOTES.md."""
    per_node = {}
    with open(ORIGINAL_N) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            doid, rid, n = int(parts[0]), int(parts[1]), float(parts[2])
            per_node.setdefault(doid, []).append((rid, n))
    return per_node


def _run_ddsim(tmp_path, doid_list, num_workers):
    """-parpath is where output gets written, and it's an absolute path --
    independent of `cwd` -- so it must point into a fresh tmp_path directory
    (never at the real, external SIPS_data tree) with just the small
    .par/.rnd/.map inputs copied in. Mirrors the instructions given for
    running this by hand; see docs/PORTING_NOTES.md."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    for ext in ("par", "rnd", "map"):
        shutil.copy(os.path.join(PARTICLES_DIR, "sips3002." + ext),
                    str(tmp_path / ("SIPS3002." + ext)))

    proc = subprocess.run(
        [sys.executable, "-m", "ddsim.DDSim", "-base", "SIPS3002",
         "-conpath", SIPS_DIR + os.sep, "-parpath", str(tmp_path) + os.sep,
         "-doid_list", ",".join(str(d) for d in doid_list), "-sv",
         "-j", str(num_workers)],
        cwd=str(tmp_path), capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-3000:]
    return n_to_exodus.node_means(n_to_exodus.read_n_file(str(tmp_path / "SIPS3002.N")))


def test_representative_subset_matches_original_via_j(tmp_path):
    orig_raw = _read_original_n()
    sentinel_doids = [d for d, rows in orig_raw.items()
                      if len(rows) == 1 and rows[0][0] == -1][:20]
    assert len(sentinel_doids) == 20, "expected >=20 sentinel doids in the real captured data"

    doids = HOT_SPOT_DOIDS + sentinel_doids
    new_means = _run_ddsim(tmp_path, doids, num_workers=2)
    orig_means = {d: sum(n for _, n in orig_raw[d]) / len(orig_raw[d]) for d in doids}

    assert set(new_means) == set(orig_means)

    # sentinel nodes: allow a rare mismatch (PORTING_NOTES.md documents ~0.01%
    # of real sentinel nodes genuinely differing, traced to the .map file on
    # disk not being a byte-exact match to whatever generated this historical
    # result -- not a code bug), but the vast majority must match exactly.
    sentinel_matches = sum(1 for d in sentinel_doids if new_means[d] == orig_means[d])
    assert sentinel_matches >= 0.9 * len(sentinel_doids), (
        f"only {sentinel_matches}/{len(sentinel_doids)} sentinel doids matched exactly")

    # hot-spot (many-particle) nodes: loose tolerance -- see PORTING_NOTES.md,
    # this is bounded by input-data provenance, not code correctness.
    for d in HOT_SPOT_DOIDS:
        reldiff = abs(new_means[d] - orig_means[d]) / orig_means[d]
        assert reldiff < 0.20, f"doid {d}: reldiff {reldiff:.3f} exceeds documented bound"


def test_j_is_deterministic_on_the_hot_spot_cluster(tmp_path):
    """The real regression guard for the multiprocessing rewrite itself
    (independent of any data-provenance question): running the same real
    doids twice, once serial once with -j, must give bit-identical results."""
    run1 = _run_ddsim(tmp_path / "run1", HOT_SPOT_DOIDS, num_workers=1)
    run2 = _run_ddsim(tmp_path / "run2", HOT_SPOT_DOIDS, num_workers=4)
    assert run1 == run2


# ---------------------------------------------------------------------------
# Whole-model (all 63,974 nodes) statistical check against a precomputed run
# ---------------------------------------------------------------------------
FULL_N = os.environ.get("DDSIM_SIPS3002_FULL_N")
full_model_skip = pytest.mark.skipif(
    not FULL_N or not os.path.exists(FULL_N),
    reason="set DDSIM_SIPS3002_FULL_N to a precomputed full-model .N file "
           "(e.g. from `ddsim -base SIPS3002 -S -sv -j <N>`) to run this -- "
           "never runs the ~1hr simulation itself")


@full_model_skip
def test_full_model_matches_original_within_documented_bounds():
    orig = n_to_exodus.node_means(n_to_exodus.read_n_file(ORIGINAL_N))
    new = n_to_exodus.node_means(n_to_exodus.read_n_file(FULL_N))

    assert set(new) == set(orig), "doid set mismatch between the two runs"
    common = sorted(orig)

    reldiffs = []
    signed = []
    for d in common:
        a, b = new[d], orig[d]
        reldiffs.append(abs(a - b) / b if b else 0.0)
        signed.append((a - b) / b if b else 0.0)

    n = len(reldiffs)
    frac_machine = sum(1 for r in reldiffs if r <= 1e-9) / n
    frac_1pct = sum(1 for r in reldiffs if r <= 0.01) / n
    frac_10pct = sum(1 for r in reldiffs if r <= 0.10) / n
    mean_signed = sum(signed) / n

    assert frac_machine >= MIN_FRAC_MACHINE_PRECISION, frac_machine
    assert frac_1pct >= MIN_FRAC_WITHIN_1PCT, frac_1pct
    assert frac_10pct >= MIN_FRAC_WITHIN_10PCT, frac_10pct
    assert abs(mean_signed) <= MAX_ABS_MEAN_SIGNED_RELDIFF, mean_signed

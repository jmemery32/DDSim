"""Golden-output regression tests against the *original* (pre-port, 2006) DDSim.

``tests/fixtures/sif_verification/`` holds a handful of small (~1000-element)
verification cubes the dissertation author built to check ``Fellipse``
(embedded) and ``Hellipse`` (surface) crack stress-intensity factors against
Abaqus, for stress fields of increasing polynomial order (constant, linear,
bilinear, biquadratic). Each ``*.output`` file is the literal captured stdout
of the original Windows/Python-2.4 DDSim.py run with ``-verify`` -- crack
lengths and K values to 6+ significant figures, plus the predicted life.

These are run through the CURRENT ``DDSim.py`` (subprocess) and every number in
the crack-growth history is compared to that captured output. This is the
strongest correctness evidence available for the port: it does not just check
that results are physically plausible, it checks that they match the actual
pre-port program, digit for digit.

One case is a known, understood exception -- see ``test_ab3333_known_divergence``.
"""
import os
import re
import shutil
import subprocess
import sys

import pytest

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "sif_verification")
DDSIM = os.path.join(os.path.dirname(__file__), "..", "DDSim.py")

# (crack_dir, mesh_basename, golden_output, doid)
CASES = [
    ("Fellipse", "constant_traction", "constant_traction.output", 666),
    ("Fellipse", "constant_traction", "constant_traction_ab3.output", 666),
    ("Fellipse", "BiLinear_traction", "BiLinear_traction.output", 666),
    ("Fellipse", "BiLinear_traction", "Bilinear_traction_ab3.output", 666),
    ("Fellipse", "Biquadratic_traction", "Biquadratic_traction.output", 666),
    ("Fellipse", "Biquadratic_traction", "biquadratic_traction_ab3.output", 666),
    ("Hellipse", "constant_traction", "constant_traction_ab1.output", 611),
    ("Hellipse", "constant_traction", "constant_traction_ab3.output", 611),
    ("Hellipse", "Linear_traction", "linear_traction_ab1.output", 611),
    ("Hellipse", "Linear_traction", "linear_traction_ab3.output", 611),
]
# see test_ab3333_known_divergence: intentionally excluded from CASES above
AB3333_CASE = ("Hellipse", "constant_traction", "constant_traction_ab3333.output", 611)

HEADER_RE = re.compile(r"Monte = (\d), a = ([\d.eE+-]+), b = ([\d.eE+-]+)")
LIFE_RE = re.compile(r"doid:\s+\d+\s+Life is:\s+([\d.eE+-]+)\s+WillGrow:\s+(-?\d+)")
BLOCK_RE = re.compile(
    r"^ (a|b|na|nb)\s+=\s+crack lengths \((\d+)\):\s*\n"
    r"([^\n]*)\n\n"
    r" Stress Intenstiy Factors \(\d+\):\s*\n"
    r"([^\n]*)\n", re.MULTILINE)
REL_TOL = 1e-4  # the golden text prints ~6-7 significant figures


def _floats(s):
    return [float(x) for x in s.replace(",", " ").split()]


def parse_run(text):
    """Extract (a, b, life, willgrow, {front: (lengths, Ks)}) from DDSim -verify output."""
    header = HEADER_RE.search(text)
    life_m = LIFE_RE.search(text)
    assert header and life_m, "could not find the expected header/life lines in:\n" + text[:2000]
    blocks = {m.group(1): (_floats(m.group(3)), _floats(m.group(4))) for m in BLOCK_RE.finditer(text)}
    return {
        "a": float(header.group(2)), "b": float(header.group(3)),
        "life": float(life_m.group(1)), "willgrow": int(life_m.group(2)),
        "blocks": blocks,
    }


def run_case(tmp_path, crack_dir, base, golden_name, doid):
    """Copy the fixture mesh, patch monte/a_b to match the golden case exactly, run DDSim.py."""
    golden = parse_run(open(os.path.join(FIXTURES, crack_dir, golden_name)).read())

    work = tmp_path / crack_dir / golden_name
    work.mkdir(parents=True)
    for ext in (".con", ".edg", ".nod", ".sig", ".smp", ".par"):
        src = os.path.join(FIXTURES, crack_dir, base + ext)
        if os.path.exists(src):
            shutil.copy(src, work / (base + ext))

    # These .par files were hand-edited across a series of exploratory runs (e.g.
    # Linear_traction.par currently holds "monte 1" with an unrelated a_b, left
    # over from a later Weibull experiment reusing the same job name). Every
    # golden .output header says "Monte = 0" with an explicit (a, b): force that,
    # rather than trusting whatever the .par currently says.
    par_path = work / (base + ".par")
    lines = par_path.read_text().splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.strip() == "monte":
            lines[i + 1] = "0\n"
        elif line.strip() == "a_b":
            lines[i + 1] = "%.10g %.10g\n" % (golden["a"], golden["b"])
    par_path.write_text("".join(lines))

    proc = subprocess.run(
        [sys.executable, DDSIM, "-base", base, "-conpath", "./", "-parpath", "./",
         "-v", "-verify", "-doid_list", str(doid), "-pa", "-sp"],
        cwd=work, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr[-3000:]
    return golden, parse_run(proc.stdout)


def assert_matches_golden(golden, new):
    assert new["life"] == pytest.approx(golden["life"], rel=REL_TOL)
    assert new["willgrow"] == golden["willgrow"]
    assert set(new["blocks"]) >= set(golden["blocks"])
    for name, (glens, gks) in golden["blocks"].items():
        nlens, nks = new["blocks"][name]
        assert nlens == pytest.approx(glens, rel=REL_TOL), "block %r crack lengths" % name
        assert nks == pytest.approx(gks, rel=REL_TOL), "block %r stress intensity factors" % name


@pytest.mark.parametrize("crack_dir,base,golden_name,doid", CASES,
                         ids=[c[2] for c in CASES])
def test_matches_original_2006_output(tmp_path, crack_dir, base, golden_name, doid):
    golden, new = run_case(tmp_path, crack_dir, base, golden_name, doid)
    assert_matches_golden(golden, new)


def test_ab3333_known_divergence(tmp_path):
    """Hellipse, a < b (Newman-Raju "a" = crack depth > "c" = surface half-length).

    The original 2006 output gives EXACTLY the same life (2035.22370851) for this
    case as for the "ab3" case (a=1.0, b=0.3333) -- to 9 significant figures.
    That would require the code to be insensitive to which of the two, very
    differently shaped, cracks was actually simulated: one crack is wide along
    the free surface and shallow (ab3), the other is narrow along the surface
    and deep (ab3333). Newman-Raju's own formula is NOT symmetric under
    depth<->surface-length exchange (only continuous at the depth == surface-length
    boundary, verified separately), and the geometric Fellipse->Hellipse
    conversion (DamClass.Fellipse._Fellipse__FtoH) does not produce swapped
    dimensions for these two inputs either -- so there is no known reason these
    two cases should agree. All available evidence points to this being a
    pre-existing bug in the ORIGINAL 2006 code (very likely a place that
    silently ignored which of a/b was larger), not something this port should
    reproduce. See PORTING_NOTES.md.

    This test locks in the port's own (internally consistent, geometry-driven)
    answer instead of the golden file, so a future change to this code path
    doesn't silently drift without anyone noticing.
    """
    crack_dir, base, golden_name, doid = AB3333_CASE
    golden, new = run_case(tmp_path, crack_dir, base, golden_name, doid)
    assert new["life"] != pytest.approx(golden["life"], rel=REL_TOL), (
        "the ab3333 divergence from the original seems to have disappeared -- "
        "if this is now expected, promote this case into CASES above")
    assert new["life"] == pytest.approx(2705.3649690778784, rel=1e-9)
    assert new["willgrow"] == 0

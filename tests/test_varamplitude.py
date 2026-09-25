import pytest

from ddsim.VarAmplitude import Spectrum


def write(tmp_path, lines):
    p = tmp_path / "spec.val"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_pairs_are_sorted_numerically(tmp_path):
    # The 2007 code sorted the *strings*, so "10.0" < "9.0" flipped (min,max).
    spec = Spectrum(write(tmp_path, ["10.0 9.0", "0.5 -0.25"]))
    assert spec.load == [(9.0, 10.0), (-0.25, 0.5)]


def test_constant_amplitude_delta_and_R(tmp_path):
    spec = Spectrum(write(tmp_path, ["0.1 1.0"] * 5))
    spec.Initialize()
    for _ in range(12):  # wraps around the end of the list
        dk, R = spec.Delta()
        assert dk == pytest.approx(0.9)
        assert R == pytest.approx(0.1)


def test_zigzag_is_condensed(tmp_path):
    # Kmax,i == Kmin,i+1  ->  the two cycles merge into a single one
    spec = Spectrum(write(tmp_path, ["0.0 0.5", "0.5 1.0", "0.0 0.5", "0.0 0.5"]))
    spec.Initialize()
    dk, R = spec.Delta()
    assert dk == pytest.approx(1.0)
    assert R == pytest.approx(0.0)
    assert spec.length == 3

"""NASGRO / Willenborg crack growth (MydadN.dadN), 7075-T651 from Table 1 of
Emery et al., Eng. Fract. Mech. 76 (2009).  Units: MPa, mm.

Table 1 gives Kie, Kic, DK1, C, n, p, q, Cthm, Cthp, Ak, Bk, Smax/sigma0,
alpha, a0, n', damp, Rso.  UTS, YS, Rcl and R are NOT in the paper's table:
the values below are typical 7075-T651 numbers and are ASSUMED.
"""
import math

import pytest

from MydadN import dadN, Willenborg

UTS, YS, RCL, R = 572.0, 503.0, 0.7, 0.1  # assumed
MATERIAL = [UTS, YS, 1320.0, 972.7, 0.75, 2.0, 3.0e-8, 2.8, 0.5, 1.0, 24.32,
            1.3, 0.1, RCL, 2.0, 0.3, R, 0.00381, 3.5, 1000.0, 3.0]


@pytest.fixture
def fat():
    return dadN(MATERIAL, 0)


def test_newman_A0_hand_calc(fat):
    # A0 = (0.825 - 0.34a + 0.05a^2) * cos(pi/2 * Smax/sigma0)^(1/a), a = 2
    expected = (0.825 - 0.68 + 0.2) * math.cos(math.pi / 2 * 0.3) ** 0.5
    assert fat.A0 == pytest.approx(expected)
    assert fat.A0 == pytest.approx(0.3257, abs=1e-4)


def test_crack_opening_function(fat):
    assert fat.Setf(0.0) == pytest.approx(fat.A0)           # R = 0
    assert fat.Setf(-0.5) == pytest.approx(fat.A0 + fat.A1 * -0.5)
    assert fat.Setf(1.0) == pytest.approx(1.0)              # f(1) = A0+A1+A2+A3 = 1


def test_Reff_small_and_large_N(fat):
    assert fat.CalcReff(0, R) == pytest.approx(0.7)         # brand-new crack
    assert fat.CalcReff(1e9, R) == pytest.approx(R)         # N >> damp
    assert fat.CalcReff(10, 0.85) == 0.85                   # R > 0.7 untouched


def test_zero_below_threshold(fat):
    a, N = 0.05, 1e6
    dkth = fat.Calc_DKth(a, fat.CalcReff(N, R))
    assert fat.Calc_dadN(0.5 * dkth, a, N, R) == 0.0


def test_growth_rate_monotonic_in_DK(fat):
    a, N = 0.05, 1e6
    dkth = fat.Calc_DKth(a, fat.CalcReff(N, R))
    rates = [fat.Calc_dadN(dkth * s, a, N, R) for s in (1.5, 2, 4, 8, 16)]
    assert all(r > 0 for r in rates)
    assert rates == sorted(rates)


def test_paris_slope_in_midrange(fat):
    # Well away from threshold and fracture the log-log slope tends to n = 2.8.
    a, N = 0.05, 1e6
    d1, d2 = 100.0, 120.0
    slope = (math.log(fat.Calc_dadN(d2, a, N, R)) - math.log(fat.Calc_dadN(d1, a, N, R))) \
        / (math.log(d2) - math.log(d1))
    assert slope == pytest.approx(fat.n, rel=0.15)


def test_threshold_lower_for_short_cracks(fat):
    # small-crack modification (NASGRO eqn 2.12): DK1* grows with a
    assert fat.Calc_DKth(0.001, 0.1) < fat.Calc_DKth(1.0, 0.1)


def test_willenborg_overload_retards_growth():
    fat = Willenborg(MATERIAL, 0)
    a, N = 0.5, 1e6
    DK = 80.0
    base = fat.Calc_dadN(DK, a, N, R)
    fat.SetKol()
    fat.Willenborg(DK, a, N, R)                       # establish baseline Kol
    fat.Willenborg(2.0 * DK, a, N, R)                 # overload cycle
    retarded = fat.Willenborg(DK, a + 1e-4, N, R)     # just after, inside plastic zone
    assert 0.0 <= retarded < base


def test_willenborg_no_overload_is_plain_nasgro():
    fat = Willenborg(MATERIAL, 0)
    a, N, DK = 0.5, 1e6, 80.0
    fat.SetKol()
    assert fat.Willenborg(DK, a, N, R) == pytest.approx(fat.Calc_dadN(DK, a, N, R))


def test_dDKth_dReff_matches_finite_difference():
    fat = Willenborg(MATERIAL, 0)
    a, Reff, h = 0.5, 0.3, 1e-6
    fd = (fat.Calc_DKth(a, Reff + h) - fat.Calc_DKth(a, Reff - h)) / (2 * h)
    assert fat.dDKthbydReff(a, Reff) == pytest.approx(fd, rel=1e-4)

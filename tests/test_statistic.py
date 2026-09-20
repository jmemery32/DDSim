import pytest

import Statistic

LIST1 = [1.2, 3.0, 3.1, 3.4, 5.4]
LIST2 = [1253.4, 1485.1, 6984.2, 4512.2, 7845.3]


def test_stata_matches_documented_values():
    stat = Statistic.Stata(2, 5)
    stat.UpdateSet()
    for a in LIST1:
        stat.UpdateSamples(a, 0)
    stat.UpdateSet()
    for a in LIST2:
        stat.UpdateSamples(a, 1)

    # reference values quoted in the original TestStata()
    assert stat.SampleMean(0) == pytest.approx(3.22)
    assert stat.SampleVariance(0) == pytest.approx(2.232)
    assert stat.SampleMean(1) == pytest.approx(4416.04)
    assert stat.SampleVariance(1) == pytest.approx(9239304.43, rel=1e-8)
    assert stat.CalcProb(0, 3.2) == pytest.approx(0.6)
    assert stat.CalcProb(1, 3000.0) == pytest.approx(0.4)


def test_stata_rid_mapping_and_duplicates():
    stat = Statistic.Stata(1, 3)
    stat.UpdateSet()
    stat.UpdateSamples(1.0, 0, 5)
    stat.UpdateSamples(1.0, 0, 7)  # same a_i -> both RIDs kept under one key
    assert stat.sample(0)[1.0] == [5, 7]
    assert stat.GetRID(1.0) == [5, 7]
    assert stat.GetRID(9.9) == -12


def test_statn_dynamic_sample_count():
    stat = Statistic.StatN(1, 'dynamic')
    stat.UpdateSet()
    for rid, n in enumerate([10.0, 20.0, 30.0], start=1):
        stat.UpdateSamples(n, 0, rid)
    assert stat.samples == 3
    assert stat.SampleMean(0) == pytest.approx(20.0)


def test_repr_is_sorted():
    stat = Statistic.Stata(1, 2)
    stat.UpdateSet()
    stat.UpdateSamples(5.0, 0)
    stat.UpdateSamples(1.0, 0)
    text = repr(stat)
    assert text.index("1.0000e+00") < text.index("5.0000e+00")

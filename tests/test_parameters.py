import pytest

from ddsim.Parameters import Parameters

PAR = """# a comment
a_b
1.0e-3 1.0e-3 2.0e-3 2.0e-3
material
{material}
dN
100.0
max_error
0.01
N_max
100000
min_inc
1.0e-6
samples
7
sets
3
monte
1
ratio
1.25
MaxCrack
5.0
"""


def test_reads_par_file(tmp_path):
    material = " ".join(str(float(i)) for i in range(21))
    (tmp_path / "job.par").write_text(PAR.format(material=material))
    p = Parameters("job", str(tmp_path) + "/")

    assert p.a_b == [[1.0e-3, 1.0e-3], [2.0e-3, 2.0e-3]]
    assert p.material == [float(i) for i in range(21)]
    assert (p.dN, p.max_error, p.N_max, p.min_inc) == (100.0, 0.01, 100000.0, 1.0e-6)
    assert (p.samples, p.sets, p.monte) == (7, 3, 1)
    assert p.r == 1.25
    assert p.Max_crack_size == 5.0
    assert "UTS" in repr(p)


def test_ratio_below_one_falls_back(tmp_path):
    material = " ".join("1.0" for _ in range(21))
    text = PAR.format(material=material).replace("ratio\n1.25", "ratio\n0.5")
    (tmp_path / "job.par").write_text(text)
    assert Parameters("job", str(tmp_path) + "/").r == 1.1

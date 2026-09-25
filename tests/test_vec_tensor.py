import math

import numpy as np
import pytest

from ddsim import JohnsVectorTools as JVT
from ddsim import Vec3D
from ddsim.ColTensor import ColTensor
from ddsim.Vec3D import CrossProd, TripleProd


def test_vec3d_arithmetic():
    a, b = Vec3D.Vec3D(1, 2, 3), Vec3D.Vec3D(4, 5, 6)
    assert list(a + b) == [5, 7, 9]
    assert list(b - a) == [3, 3, 3]
    assert a * b == 32.0                       # dot product, a float
    assert list(2 * a) == [2, 4, 6]
    assert list(a * 2.5) == [2.5, 5.0, 7.5]
    assert list(a / 2) == [0.5, 1.0, 1.5]
    assert list(-a) == [-1, -2, -3]
    assert list(a + 1) == [2, 3, 4]            # scalar applies to every component
    assert list(10 - a) == [9, 8, 7]
    with pytest.raises(TypeError):
        a * "x"


def test_vec3d_methods():
    v = Vec3D.Vec3D(3, 4, 0)
    assert v.Magnitude() == 5.0
    n = v.Normalize()
    assert (n.x(), n.y(), n.z()) == pytest.approx((0.6, 0.8, 0.0))
    assert list(v) == [3, 4, 0]                # Normalize returns a NEW vector
    with pytest.raises(TypeError):
        Vec3D.Vec3D(0, 0, 0).Normalize()
    assert repr(Vec3D.Vec3D(1, 2.5, -3)) == "1 2.5 -3"
    m = Vec3D.Vec3D(1, 5, 3)
    m.Max(Vec3D.Vec3D(2, 2, 2))
    assert list(m) == [2, 5, 3]
    m.Min(Vec3D.Vec3D(0, 9, 1))
    assert list(m) == [0, 5, 1]


def test_vec3d_inplace_mutates():
    v = Vec3D.Vec3D(1, 1, 1)
    alias = v
    v += Vec3D.Vec3D(1, 2, 3)
    v *= 2
    assert alias is v and list(v) == [4, 6, 8]


def test_cross_and_triple():
    x, y, z = Vec3D.Vec3D(1, 0, 0), Vec3D.Vec3D(0, 1, 0), Vec3D.Vec3D(0, 0, 1)
    assert list(CrossProd(x, y)) == [0, 0, 1]
    assert TripleProd(x, y, z) == 1.0


def test_coltensor_components_and_invariants():
    t = ColTensor(1, 2, 3, 0.4, 0.5, 0.6)      # xx yy zz xy yz zx
    assert (t.xx(), t.yy(), t.zz(), t.xy(), t.yz(), t.zx()) == (1, 2, 3, 0.4, 0.5, 0.6)
    assert t.I1() == 6
    assert t.I3() == pytest.approx(np.linalg.det(t.matrix()))
    uni = ColTensor(100, 0, 0, 0, 0, 0)
    assert uni.EffectiveStress() == pytest.approx(100.0)     # uniaxial: von Mises = sigma


def test_principal_values_sorted_with_eigenvectors():
    t = ColTensor(10, 30, 20, 0, 0, 0)
    prince, evect = t.PrincipalValues()
    assert (prince.xx(), prince.yy(), prince.zz()) == pytest.approx((30, 20, 10))
    e1 = evect[0]
    assert abs(e1.y()) == pytest.approx(1.0)                 # sigma1 = syy -> y axis

    rng = np.random.default_rng(1)
    for _ in range(20):
        t = ColTensor(*rng.normal(size=6) * 100)
        prince, evect = t.PrincipalValues()
        s = [prince.xx(), prince.yy(), prince.zz()]
        assert s[0] >= s[1] >= s[2]
        for lam, e in zip(s, evect):
            assert t.matrix() @ np.array(list(e)) == pytest.approx(lam * np.array(list(e)))
            assert e.Magnitude() == pytest.approx(1.0)


def test_jvt():
    assert JVT.Plus([1, 2], [3, 4]) == [4, 6]
    assert JVT.Minus([5, 5], [1, 2]) == [4, 3]
    assert JVT.Star([2, 3], [4, 5]) == [8, 15]
    assert JVT.Divide([8, 9], [2, 3]) == [4, 3]
    assert JVT.ScalarMult([1, 2], 3) == [3, 6]


def test_numpy_scalars_defer_to_vec3d_operators():
    v = Vec3D.Vec3D(1, 2, 3)
    for k in (np.float64(2.0), np.float32(2.0), np.int64(2)):
        for r in (k * v, v * k):
            assert isinstance(r, Vec3D.Vec3D) and list(r) == [2, 4, 6]
    assert isinstance(np.float64(1.0) + v, Vec3D.Vec3D)
    t = ColTensor(1, 2, 3, 4, 5, 6)
    assert isinstance(np.float64(2.0) * t, ColTensor)

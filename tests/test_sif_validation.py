"""End-to-end validation of the stress intensity factors against closed-form solutions.

example1 is a 20 x 20 x 20 cube with a uniform stress field (syy = 3.4), so
textbook crack solutions apply exactly:

* an *embedded* elliptical crack (Irwin):       K = s sqrt(pi b) / E(k) * (sin^2 + (b/a)^2 cos^2)^(1/4)
* a *surface* semi-circular crack (Newman-Raju, a/t -> 0)

These exercise mesh interpolation, principal-stress rotation, Vec3D/ColTensor,
GeomUtils (surface detection) and the SIF machinery together.
"""
import math
import os

import numpy as np
import pytest
from scipy.special import ellipe

import DamClass
import DamMo
import MeshTools
import Parameters

EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "examples", "example1") + os.sep
SIGMA = 3.4


@pytest.fixture(scope="module")
def setup():
    model = MeshTools.MeshTools(EXAMPLE + "example1", "RDB")
    model.SetPointInsideTolerance(1.0e-7)
    params = Parameters.Parameters("example1", EXAMPLE)
    cracks = DamMo.DamModel(model, model.GetNodeList(), 0, 0, "0", 0, params, None, None, None)
    return model, params, cracks


def crack_at(setup, nid, a, b):
    model, params, cracks = setup
    xyz = model.GetNodeInfo(nid)[0]
    fe = DamClass.Fellipse(nid, xyz, model, a, b, params.material, None, 0, 0)
    changed, _ = fe.GeometryCheck(cracks.CMesh)
    return fe if changed == 0 else changed


@pytest.mark.parametrize("a", [0.25, 0.5, 2.0])
def test_embedded_penny_crack(setup, a):
    K = crack_at(setup, 8, a, a).CalculateKi(1.0)
    assert list(K) == pytest.approx([2 * SIGMA * math.sqrt(a / math.pi)] * 4, rel=1e-4)


@pytest.mark.parametrize("a,b", [(1.0, 0.5), (0.5, 1.0), (2.0, 1.0), (1.0, 4.0)])
def test_embedded_elliptical_crack_matches_irwin(setup, a, b):
    amaj, bmin = max(a, b), min(a, b)
    Ebk = ellipe(1.0 - (bmin / amaj) ** 2)
    k_minor_end = SIGMA * math.sqrt(math.pi * bmin) / Ebk
    k_major_end = k_minor_end * math.sqrt(bmin / amaj)
    el = crack_at(setup, 8, a, b)
    assert type(el).__name__ == "Fellipse"
    K = el.CalculateKi(1.0)
    # K[0], K[2] are at the two ends of the `a` axis, K[1], K[3] at the `b` axis
    a_end, b_end = (k_major_end, k_minor_end) if a >= b else (k_minor_end, k_major_end)
    assert [K[0], K[1], K[2], K[3]] == pytest.approx([a_end, b_end, a_end, b_end], rel=1e-4)


def newman_raju_semicircle(a, phi, scale=1.0):
    """Semicircular surface crack in a thick plate (a/c = 1, a/t -> 0)."""
    Q = 1.0 + 1.464
    M1 = 1.13 - 0.09
    g = 1.0 + 0.1 * (1.0 - math.sin(phi)) ** 2
    return scale * SIGMA * math.sqrt(math.pi * a / Q) * M1 * g


@pytest.mark.parametrize("a", [0.01, 0.5, 1.0])
def test_surface_crack_matches_newman_raju(setup, a):
    scale = 60.0
    el = crack_at(setup, 10, a, a)                       # face-centre node (0, 0, 10)
    assert type(el).__name__ == "Hellipse"               # detected as a surface crack
    K = el.CalculateKi(scale)
    surface = newman_raju_semicircle(a, 0.0, scale)
    deepest = newman_raju_semicircle(a, math.pi / 2, scale)
    assert list(K) == pytest.approx([surface, deepest, surface], rel=2e-3)


def test_interior_node_is_not_a_surface_crack(setup):
    assert type(crack_at(setup, 8, 0.5, 0.5)).__name__ == "Fellipse"

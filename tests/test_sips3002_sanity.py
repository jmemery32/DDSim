"""Sanity checks against the real SIPS3002 open-hole validation model
(Emery et al., "DDSim... Part I", Eng. Fract. Mech. 76 (2009) 1500-1530, Sec. 6.1).

This is the actual proprietary NGC coupon geometry from the dissertation
validation study, kept OUTSIDE the repository (see PORTING_NOTES.md). Point
``DDSIM_SIPS3002_DIR`` at a local copy of
``SIPS_data/DDSimLI/SIPS3002_open/ConstantAmplitude`` (containing
``SIPS3002.con/.edg/.nod/.sig/.smp``) to run these; otherwise they are skipped.
These do not touch crack growth at all -- just the mesh reader and surface
detection -- so they need no `.par`/`.rnd`/`.map` and run in a few seconds.
"""
import os

import numpy as np
import pytest

from ddsim import MeshTools

SIPS_DIR = os.environ.get("DDSIM_SIPS3002_DIR")
pytestmark = pytest.mark.skipif(
    not SIPS_DIR or not os.path.exists(os.path.join(SIPS_DIR, "SIPS3002.con")),
    reason="set DDSIM_SIPS3002_DIR to a local copy of the (proprietary, non-repo) "
           "SIPS3002_open/ConstantAmplitude mesh to run these")


@pytest.fixture(scope="module")
def model():
    return MeshTools.MeshTools(os.path.join(SIPS_DIR, "SIPS3002"), "RDB")


def test_element_and_node_counts_match_the_paper(model):
    """Sec. 6.1: "140,024 8-noded brick elements and 184 6-noded wedge elements
    ... 63,974 nodes on the surface of the model." """
    types = [model._eclass[i].name for i in range(len(model._eclass))]
    assert types.count("BRICK_8") == 140024
    assert types.count("WEDGE_6") == 184
    assert len(types) == 140208
    assert len(model.GetNodeList()) == 172601
    surface_nodes = sum(1 for n in model.GetNodeList() if model.IsSurfaceNode(n))
    assert surface_nodes == 63974


def test_peak_stress_forms_a_single_tight_cluster(model):
    """The paper reports one dominant hot-spot ("the aft side of hole 16 near the
    intersection of the counter-bore with the main-bore"), not several unrelated
    concentrations. Take the highest-sigma_yy nodes and check they are all within
    a small ball of each other -- a real, localized stress concentration, not
    scattered numerical noise (which would indicate a mesh-reading or
    stress-alignment bug on a model this large)."""
    nodes = model.GetNodeList()
    syy = np.array([model.GetNodeInfo(n)[2].yy() for n in nodes])
    top = np.argsort(syy)[-25:]                       # top 25 of 172,601 nodes
    coords = np.array([list(model.GetNodeInfo(nodes[i])[0]) for i in top])

    assert syy[top].min() > 0.5 * syy[top].max()       # a real peak, not one outlier node
    spread = coords.max(axis=0) - coords.min(axis=0)
    assert np.all(spread < 2.0), "top-stress nodes are not spatially clustered: %s" % spread
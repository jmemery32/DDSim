"""Mesh + nodal-stress input for DDSim Level I.

``MeshData`` is the format-neutral container the query engine (``MeshTools``)
consumes.  Readers produce it:

* :func:`read_rdb`  -- the Cornell "relational DB" ASCII files (``.smp .nod .con
  .edg .sig``) used by the 2007 code.
* (planned) an Exodus II reader producing the same ``MeshData``.

Conventions
-----------
* ``connectivity[i]`` lists the node ids of element ``i`` in the node order of
  the corresponding element in :mod:`elements` (corners first, then mid-side
  nodes), i.e. exactly what the 2007 ``MeshTools`` built internally.
* ``stress`` rows are ``(sxx, syy, szz, sxy, syz, szx)``, one per node, aligned
  with ``node_ids``.
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class MeshData:
    node_ids: np.ndarray                  # (n,) int
    coords: np.ndarray                    # (n, 3) float
    elem_ids: np.ndarray                  # (m,) int
    elem_types: list                      # (m,) str, keys of elements.ELEMENTS
    connectivity: list                    # (m,) arrays of node ids
    stress: np.ndarray = None             # (n, 6) float, aligned with node_ids
    meta: dict = field(default_factory=dict)

    def validate(self):
        n = len(self.node_ids)
        if self.coords.shape != (n, 3):
            raise ValueError("coords must be (n_nodes, 3)")
        if self.stress is not None and self.stress.shape != (n, 6):
            raise ValueError("stress must be (n_nodes, 6)")
        if len(set(self.node_ids.tolist())) != n:
            raise ValueError("duplicate node ids")
        known = set(self.node_ids.tolist())
        for eid, conn in zip(self.elem_ids, self.connectivity):
            missing = [c for c in conn if int(c) not in known]
            if missing:
                raise ValueError("element %d references unknown nodes %s" % (eid, missing[:5]))


# ---------------------------------------------------------------------------
# RDB reader
# ---------------------------------------------------------------------------
_SHAPE_TO_TYPE = {
    "Tet1IShapeFunc": "TET_4",
    "Tet2IShapeFunc": "TET_10",
    "Tet2I2FShapeFunc": "TET_10",
    "Brick1IShapeFunc": "BRICK_8",
    "Brick2IShapeFunc": "BRICK_20",
    "Wedge1IShapeFunc": "WEDGE_6",
    "Wedge2IShapeFunc": "WEDGE_15",
}
_UNSUPPORTED = {"Pyramid1IShapeFunc", "Pyramid2IShapeFunc"}

# local edge (pair of corner indices) of each mid-side node, in mid-side node order
_MIDSIDE_EDGES = {
    "TET_10": [(0, 1), (1, 2), (2, 0), (2, 3), (1, 3), (0, 3)],
    "BRICK_20": [(0, 1), (2, 3), (4, 5), (6, 7), (0, 2), (1, 3), (4, 6), (5, 7),
                 (0, 4), (1, 5), (2, 6), (3, 7)],
    "WEDGE_15": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
}


def _rows(path):
    """Whitespace-split non-blank, non-comment lines of a text file."""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                yield line.split()


def read_rdb(base):
    """Read ``base``.smp/.nod/.con/.sig (and .edg for quadratic elements)."""
    shape_type = {}
    quadratic_ids = set()
    for row in _rows(base + ".smp"):
        sid, name = int(row[0]), row[1]
        if name in _UNSUPPORTED:
            raise NotImplementedError("%s elements are not supported" % name)
        etype = _SHAPE_TO_TYPE[name]
        shape_type[sid] = etype
        if etype in _MIDSIDE_EDGES:
            quadratic_ids.add(sid)

    node_ids, coords = [], []
    for row in _rows(base + ".nod"):
        node_ids.append(int(row[0]))
        coords.append([float(row[1]), float(row[2]), float(row[3])])

    elems = []                                   # (eid, shape id, corner node ids)
    need_edges = False
    for row in _rows(base + ".con"):
        eid, sid, ncorner = int(row[0]), int(row[1]), int(row[3])
        elems.append((eid, sid, [int(v) for v in row[4:4 + ncorner]]))
        need_edges = need_edges or sid in quadratic_ids

    side_node = {}
    if need_edges:
        for row in _rows(base + ".edg"):
            a, b = int(row[3]), int(row[4])
            side_node[(min(a, b), max(a, b))] = int(row[0])

    elem_ids, elem_types, connectivity = [], [], []
    for eid, sid, corners in elems:
        etype = shape_type[sid]
        conn = list(corners)
        for a, b in _MIDSIDE_EDGES.get(etype, ()):
            na, nb = corners[a], corners[b]
            conn.append(side_node[(min(na, nb), max(na, nb))])
        elem_ids.append(eid)
        elem_types.append(etype)
        connectivity.append(np.array(conn, dtype=np.int64))

    index = {nid: i for i, nid in enumerate(node_ids)}
    stress = np.zeros((len(node_ids), 6))
    for row in _rows(base + ".sig"):
        stress[index[int(row[0])]] = [float(v) for v in row[1:7]]

    data = MeshData(node_ids=np.array(node_ids, dtype=np.int64),
                    coords=np.array(coords, dtype=float),
                    elem_ids=np.array(elem_ids, dtype=np.int64),
                    elem_types=elem_types, connectivity=connectivity,
                    stress=stress, meta={"format": "RDB", "base": base})
    data.validate()
    return data

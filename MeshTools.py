"""Finite element mesh queries for DDSim.

Pure-Python/numpy replacement for the compiled ``MeshTools.pyd`` (Python 2
extension built from ``MeshTools.cpp``).  The public interface -- method names,
return values and exceptions -- is the one DDSim's ``DamClass``/``DamMo`` were
written against::

    model = MeshTools.MeshTools(path + name, 'RDB')
    model.SetPointInsideTolerance(1.0e-7)
    stress = model.GetPtStress(Vec3D(x, y, z))       # ColTensor
    status, dist = model.IsPointOutsideMesh(Vec3D(x, y, z))

Point-query status codes (``IsPointOutsideMesh`` / internal ``_pt_query``)::

    -3   no element bounding box contains the point: clearly outside
    -2   point is inside an element
    -1   point is outside, but the nearest element/point (and distance) is known
    >=0  Newton iteration for natural coordinates did not converge (element id)

Deliberate differences from the 2007 implementation
---------------------------------------------------
* Surface normals are oriented **outward** using the parent element (the 2007
  face tables gave inward normals for right-handed 6-node wedges).
* ``NearestPoint`` for simplex elements is a correct projection (see
  :mod:`elements`).
* The range tree is replaced by a uniform grid over the element bounding boxes.
* The nodal displacement field is not stored (it was a hard-wired zero vector).
"""
import numpy as np

import ColTensor as _CT
import Vec3D as _V
import _kernels as _K
import elements as _E
import mesh_io as _io


class EmptySearchResult(Exception):
    """A point query found no usable element."""


class InvalidNodeId(Exception):
    """Unknown node id."""


class InvalidElemId(Exception):
    """Unknown element id."""


_MAX_ITS = 12          # Newton iterations for the natural coordinates
_BOX_PAD = 0.1         # element bounding boxes are grown by 10% of their diagonal


class MeshTools:

    def __init__(self, name=None, model_type='RDB', data=None):
        self.PointInsideTol = 1.0e-7
        self._cached = None            # index of the element the last query landed in
        self._max_min = None
        if data is None:
            if model_type == 'RDB':
                data = _io.read_rdb(name)
            else:
                raise ValueError("unsupported model type %r" % (model_type,))
        self._load(data)

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    def _load(self, data):
        self._node_ids = [int(n) for n in data.node_ids]
        self._row = {nid: i for i, nid in enumerate(self._node_ids)}
        self._xyz = np.asarray(data.coords, dtype=float)
        stress = data.stress
        self._sig = np.zeros((len(self._node_ids), 6)) if stress is None else np.asarray(stress, dtype=float)

        self._elem_ids = [int(e) for e in data.elem_ids]
        self._eindex = {eid: i for i, eid in enumerate(self._elem_ids)}
        self._eclass = [_E.ELEMENTS[t] for t in data.elem_types]
        self._enodes = [[int(n) for n in c] for c in data.connectivity]
        self._erows = [np.array([self._row[n] for n in c]) for c in self._enodes]

        # padded arrays for the compiled point-location kernels
        m = len(self._erows)
        self._code = np.array([e.code for e in self._eclass], dtype=np.int64)
        self._conn = np.zeros((m, _K.MAX_NODES), dtype=np.int64)
        for i, rows in enumerate(self._erows):
            self._conn[i, :len(rows)] = rows
        self._xyz = np.ascontiguousarray(self._xyz)
        self._sig = np.ascontiguousarray(self._sig)
        # scratch buffers (a MeshTools object is not thread-safe; processes are fine)
        self._nc = np.empty(3)
        self._near = np.empty(3)
        self._X = np.empty((_K.MAX_NODES, 3))
        self._N = np.empty(_K.MAX_NODES)
        self._dN = np.empty((_K.MAX_NODES, 3))
        self._tmp = np.empty(_K.MAX_NODES, dtype=np.complex128)
        self._out6 = np.empty(6)

        self._build_bounding_boxes()
        self._build_reverse_connectivity()
        self._build_surface_mesh()

    def _build_bounding_boxes(self):
        m = len(self._erows)
        lo, hi = np.empty((m, 3)), np.empty((m, 3))
        for i, rows in enumerate(self._erows):
            pts = self._xyz[rows]
            lo[i], hi[i] = pts.min(axis=0), pts.max(axis=0)
        pad = _BOX_PAD * np.linalg.norm(hi - lo, axis=1)[:, None]
        self._lo, self._hi = lo - pad, hi + pad

        # uniform grid over the boxes (compiled CSR layout, see _kernels.build_grid)
        extent = (self._hi - self._lo).max(axis=1)
        self._origin = np.ascontiguousarray(self._lo.min(axis=0))
        cell = float(np.median(extent)) if m else 1.0
        span = self._hi.max(axis=0) - self._origin
        while np.prod(np.floor(span / cell) + 1.0) > 2.0e7:      # keep the grid small
            cell *= 1.5
        self._inv_cell = 1.0 / cell
        self._dims = (np.floor(span * self._inv_cell).astype(np.int64) + 1)
        self._lo = np.ascontiguousarray(self._lo)
        self._hi = np.ascontiguousarray(self._hi)
        self._cell_start, self._cell_items = _K.build_grid(
            self._lo, self._hi, self._origin, self._inv_cell, self._dims)

    def _build_reverse_connectivity(self):
        self._rev = {}
        for i, nodes in enumerate(self._enodes):
            for n in nodes:
                self._rev.setdefault(n, []).append(self._elem_ids[i])

    def _build_surface_mesh(self):
        """Faces that belong to exactly one element form the surface mesh."""
        faces = {}
        for i, (el, nodes) in enumerate(zip(self._eclass, self._enodes)):
            for face in el.faces:
                fn = [nodes[k] for k in face]
                ncorner = len(fn) if len(fn) <= 4 else len(fn) // 2
                key = tuple(sorted(fn[:ncorner]))
                if key in faces:
                    del faces[key]
                else:
                    faces[key] = (fn, i)

        self._surf_nodes = []          # node ids of each surface element
        self._surf_parent = []         # index of the element the face belongs to
        self._surf_rev = {}            # node id -> surface element ids
        for seid, (fn, parent) in enumerate(faces.values()):
            self._surf_nodes.append(fn)
            self._surf_parent.append(parent)
            for n in fn:
                self._surf_rev.setdefault(n, []).append(seid)

    # ------------------------------------------------------------------
    # small accessors
    # ------------------------------------------------------------------
    def SetPointInsideTolerance(self, tol):
        self.PointInsideTol = float(tol)

    def GetNodeList(self):
        return list(self._node_ids)

    def GetNodeInfo(self, nid):
        """-> (coordinates Vec3D, displacement Vec3D (zero), stress ColTensor)."""
        try:
            i = self._row[int(nid)]
        except KeyError:
            raise InvalidNodeId("unknown node ID")
        x = self._xyz[i]
        return (_V.Vec3D(x[0], x[1], x[2]), _V.Vec3D(0.0, 0.0, 0.0), _CT.ColTensor(*self._sig[i]))

    def GetElemInfo(self, eid):
        try:
            return list(self._enodes[self._eindex[int(eid)]])
        except KeyError:
            raise InvalidElemId("unknown element ID")

    def GetAdjacentElems(self, nid):
        try:
            return list(self._rev[int(nid)])
        except KeyError:
            raise InvalidNodeId("unknown node ID")

    def GetSurfElemInfo(self, seid):
        try:
            return list(self._surf_nodes[int(seid)])
        except IndexError:
            raise InvalidElemId("unknown element ID")

    def GetAdjacentSurfElems(self, nid):
        try:
            return list(self._surf_rev[int(nid)])
        except KeyError:
            raise InvalidNodeId("unknown node ID")

    def IsSurfaceNode(self, nid):
        return 1 if int(nid) in self._surf_rev else 0

    def GetMaxDimension(self):
        """-> (Max, Min): opposite corners of the model's bounding box (Vec3D)."""
        if self._max_min is None:
            mx, mn = self._xyz.max(axis=0), self._xyz.min(axis=0)
            self._max_min = (_V.Vec3D(*mx), _V.Vec3D(*mn))
        return self._max_min

    # ------------------------------------------------------------------
    # surface normals
    # ------------------------------------------------------------------
    def _one_surface_normal(self, nid, seid):
        fn = self._surf_nodes[seid]
        local = fn.index(nid)
        coords = self._xyz[[self._row[n] for n in fn]]
        j0, j1 = _E.surface_node_tangents(len(fn), local, coords)
        n = np.cross(j0, j1)
        # orient away from the element the face belongs to
        parent = self._xyz[self._erows[self._surf_parent[seid]]].mean(axis=0)
        if n @ (coords.mean(axis=0) - parent) < 0.0:
            n = -n
        return n / np.linalg.norm(n)

    def SurfaceNormal(self, nid, eid=-1):
        """Outward unit normal at surface node ``nid``.

        With ``eid`` the normal of that one surface element, otherwise the
        (re-normalised) average over the adjacent surface elements.
        """
        nid = int(nid)
        if eid is not None and eid >= 0:
            n = self._one_surface_normal(nid, int(eid))
        else:
            try:
                seids = self._surf_rev[nid]
            except KeyError:
                raise InvalidNodeId("unknown node ID")
            total = sum(self._one_surface_normal(nid, s) for s in seids)
            n = total / np.linalg.norm(total)
        return _V.Vec3D(n[0], n[1], n[2])

    # ------------------------------------------------------------------
    # point queries
    # ------------------------------------------------------------------
    @staticmethod
    def _as_point(pt):
        return np.array([pt[0], pt[1], pt[2]], dtype=float)

    def _pt_query(self, q):
        """-> (status, element index, natural coords, distance); see the module docstring."""
        cached = -1 if self._cached is None else self._cached
        nat = np.empty(3)
        status, i, dist = _K.pt_query(
            q, cached, self.PointInsideTol, _MAX_ITS, self._origin, self._inv_cell,
            self._dims, self._cell_start, self._cell_items, self._lo, self._hi,
            self._code, self._conn, self._xyz, nat, self._nc, self._X, self._N,
            self._dN, self._tmp, self._near)
        if status == -2:
            self._cached = i
        elif status == -1:
            self._cached = None
        return status, i, nat, dist

    def IsPointOutsideMesh(self, pt):
        """-> (status, distance); see the module docstring for the status codes."""
        status, i, _, dist = self._pt_query(self._as_point(pt))
        if status == -4:                    # Newton did not converge: report the element id
            return self._elem_ids[i], 0.0
        return status, dist

    def GetPtStresses(self, points):
        """Stress at many points at once: ``(n, 3)`` array-like -> ``(n, 6)`` ndarray
        with columns (sxx, syy, szz, sxy, syz, szx).  Raises ``EmptySearchResult`` if
        any point is not in, or just next to, the mesh."""
        pts = np.ascontiguousarray(points, dtype=float).reshape(-1, 3)
        out = np.empty((len(pts), 6))
        cached = -1 if self._cached is None else self._cached
        done, cached = _K.stress_batch(
            pts, cached, self.PointInsideTol, _MAX_ITS, self._origin, self._inv_cell,
            self._dims, self._cell_start, self._cell_items, self._lo, self._hi,
            self._code, self._conn, self._xyz, self._sig, out)
        self._cached = None if cached < 0 else int(cached)
        if done < len(pts):
            raise EmptySearchResult("empty range tree search")
        return out

    def GetPtStress(self, pt):
        """Stress tensor (ColTensor) interpolated at ``pt``.

        Raises ``EmptySearchResult`` if the point is not in, or just next to,
        the mesh.
        """
        return _CT.ColTensor(*self.GetPtStresses(self._as_point(pt))[0])

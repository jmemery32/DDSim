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
import itertools

import numpy as np

import ColTensor as _CT
import Vec3D as _V
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
        self._cached = None            # (element index, element class, node coords)
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

        # uniform grid: each element is registered in every cell its box touches
        extent = (self._hi - self._lo).max(axis=1)
        self._origin = self._lo.min(axis=0)
        self._cell = float(np.median(extent)) if m else 1.0
        self._grid = {}
        ilo = np.floor((self._lo - self._origin) / self._cell).astype(int)
        ihi = np.floor((self._hi - self._origin) / self._cell).astype(int)
        for i in range(m):
            for c in itertools.product(*(range(ilo[i, k], ihi[i, k] + 1) for k in range(3))):
                self._grid.setdefault(c, []).append(i)

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

    def _candidates(self, q):
        cell = tuple(np.floor((q - self._origin) / self._cell).astype(int))
        idx = self._grid.get(cell)
        if not idx:
            return []
        idx = np.array(idx)
        inside = np.all((self._lo[idx] <= q) & (q <= self._hi[idx]), axis=1)
        return idx[inside].tolist()

    def _find_natural(self, el, X, q):
        """Newton iteration for the natural coordinates of ``q`` in an element.

        Returns the natural coordinates, or ``None`` if it did not converge.
        """
        tol = self.PointInsideTol
        nc = np.array(el.center, dtype=float)
        for _ in range(_MAX_ITS):
            guess = el.shape(nc) @ X
            jac = el.dshape(nc).T @ X            # jac[i, j] = d x_j / d nc_i
            try:
                update = np.linalg.solve(jac.T, guess - q)
            except np.linalg.LinAlgError:
                return None
            new = nc - update
            if not np.all(np.isfinite(new)):
                return None
            if np.any(np.abs(new - nc) > tol):
                nc = new
            else:
                return new
        return None

    def _pt_query(self, q):
        """-> (status, element index, element class, natural coords, distance)."""
        tol = self.PointInsideTol

        # first try the element the previous query landed in, with a negative
        # tolerance so a point close to a shared face is re-searched properly
        if self._cached is not None:
            i, el, X = self._cached
            nat = self._find_natural(el, X, q)
            if nat is not None and el.point_inside(nat, -tol):
                return -2, i, el, nat, 0.0

        cand = self._candidates(q)
        if not cand:
            return -3, None, None, None, 0.0

        rvalue = -1
        best = None                                   # (dist, element idx, class, nat)
        for i in cand:
            el = self._eclass[i]
            X = self._xyz[self._erows[i]]
            nat = self._find_natural(el, X, q)
            if nat is None:                           # did not converge
                rvalue = self._elem_ids[i]
                continue
            if el.point_inside(nat, tol):
                self._cached = (i, el, X)
                return -2, i, el, nat, 0.0
            surf = el.nearest_point(nat)
            dist = float(np.linalg.norm(q - el.shape(surf) @ X))
            if best is None or dist < best[0]:
                best = (dist, i, el, surf)

        if best is None:                              # nothing converged
            return rvalue, None, None, None, 0.0

        # the point is (numerically) just outside: report the nearest element
        self._cached = None
        dist, i, el, surf = best
        return -1, i, el, surf, dist

    def IsPointOutsideMesh(self, pt):
        """-> (status, distance); see the module docstring for the status codes."""
        status, _, _, _, dist = self._pt_query(self._as_point(pt))
        return status, dist

    def GetPtStress(self, pt):
        """Stress tensor (ColTensor) interpolated at ``pt``.

        Raises ``EmptySearchResult`` if the point is not in, or just next to,
        the mesh.
        """
        status, i, el, nat, _ = self._pt_query(self._as_point(pt))
        if status not in (-1, -2):
            raise EmptySearchResult("empty range tree search")
        return _CT.ColTensor(*(el.shape(nat) @ self._sig[self._erows[i]]))

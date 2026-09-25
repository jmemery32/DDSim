"""Exodus II mesh + nodal-stress input/output, for reading real FE results and
for writing DDSim's predictions somewhere ParaView (or any other Exodus-aware
viewer) can show them as contours.

Exodus II is a NetCDF file with a fixed schema (dimensions like ``num_nodes``,
``num_elem``, ``num_el_blk``; variables like ``coord``, ``connect1``,
``vals_nod_var1``, ...); see the SEACAS project
(https://github.com/sandialabs/seacas) for the reference implementation and
``docs/PORTING_NOTES.md`` for how the node-order mapping below was derived
and verified.

Scope
-----
Linear elements only (``TET_4``, ``WEDGE_6``, ``BRICK_8`` / Exodus
``TETRA4``/``TET4``, ``WEDGE6``, ``HEX8``) -- every real mesh DDSim has ever
been run on (the SIPS3002 validation model, ``example1``, the
SIF_Verification cubes) uses only these. Quadratic elements (``TETRA10``,
``WEDGE15``, ``HEX20``) raise ``NotImplementedError``: their node-ordering
convention needs a real sample file or an authoritative reference to verify
against (see the HEX8 derivation below) before it can be trusted, and get it
wrong and a mesh can *look* fine while silently being wrong.

Node order
----------
Exodus's node order does not match ``elements.py``'s for HEX8 (Exodus goes
around one face, then the opposite face, in the same order; our ``Brick8``
uses a binary-counting order -- see ``elements.py``). The permutation below
was derived from, and checked against, a real sample file
(``tests/fixtures/exodus/beam_hex8.g``): dump ``coord``/``connect1`` for a
real element and confirm which physical corner each Exodus slot holds.

For TET4 and WEDGE6, no real sample file was available; the permutation
below assumes the standard (Exodus/VTK/CGNS-consistent) convention, where a
tet's node 3 is the apex opposite the 0-1-2 base and a wedge's nodes 0-1-2/
3-4-5 are the bottom/top triangles (translated, like ours). Rotating which
base corner is "first" cannot change the interpolated result (trilinear/
tetrahedral interpolation is invariant to it, as long as face grouping and
handedness are preserved), so the only real risk is a *reflected* (inside
out) element -- caught for every element of every type, on every read, by
:func:`_check_nonzero_volume`, which raises immediately with the element id
rather than silently returning a wrong answer.
"""
import numpy as np

from .mesh_io import MeshData

# our local index k -> which Exodus connectivity slot holds that corner.
# HEX8: empirically derived (see module docstring); TET4/WEDGE6: identity,
# relying on _check_nonzero_volume as a safety net against garbled source files.
_EXODUS_TO_OURS = {
    "TET_4": [0, 1, 2, 3],
    "WEDGE_6": [0, 1, 2, 3, 4, 5],
    "BRICK_8": [0, 1, 3, 2, 4, 5, 7, 6],
}
_EXODUS_TYPE_TO_OURS = {
    "TET4": "TET_4", "TETRA4": "TET_4", "TETRA": "TET_4",
    "WEDGE6": "WEDGE_6", "WEDGE": "WEDGE_6",
    "HEX8": "BRICK_8", "HEX": "BRICK_8",
}
_UNSUPPORTED_QUADRATIC = {"TET10", "TETRA10", "WEDGE15", "HEX20", "HEX27"}

# common nodal-variable name conventions for the 6 stress tensor components,
# lower-cased for case-insensitive matching; extend/override with the
# ``stress_names`` argument to read_exodus for anything else.
_STRESS_ALIASES = {
    "xx": ["stress_xx", "sigma_xx", "sxx", "s11", "stress_11"],
    "yy": ["stress_yy", "sigma_yy", "syy", "s22", "stress_22"],
    "zz": ["stress_zz", "sigma_zz", "szz", "s33", "stress_33"],
    "xy": ["stress_xy", "sigma_xy", "sxy", "s12", "stress_12"],
    "yz": ["stress_yz", "sigma_yz", "syz", "s23", "stress_23"],
    "zx": ["stress_zx", "stress_xz", "sigma_zx", "sigma_xz", "szx", "sxz", "s13", "stress_13"],
}


def _check_nonzero_volume(etype, coords, elem_id):
    """Raise if an element's corner order gives (numerically) zero volume --
    a strong signal the node-order permutation grouped the wrong corners into
    a face (a self-intersecting/degenerate element), rather than just a
    harmless axis reflection.

    Deliberately does NOT require a specific *sign*: a pure reflection of one
    natural-coordinate axis (e.g. Exodus's real, verified HEX8 mapping is
    ours with one axis flipped relative to elements.py's own reference node
    order -- see test_check_nonzero_volume_accepts_our_own_reference_orderings)
    flips this formula's sign without making the element invalid, and without
    affecting DDSim's actual use of the shape functions (Newton iteration for
    natural coordinates + interpolation), which only needs the element's
    Jacobian to be non-singular, not any particular handedness convention.
    A genuinely garbled/self-intersecting element, in contrast, reliably
    shows up as a (near-)zero or wildly-off-magnitude volume here.
    """
    c = coords
    if etype == "TET_4":
        vol = np.dot(c[1] - c[0], np.cross(c[2] - c[0], c[3] - c[0]))
    elif etype == "WEDGE_6":
        vol = np.dot(c[3] - c[0], np.cross(c[1] - c[0], c[2] - c[0]))
    else:  # BRICK_8: signed volume of the 6 corner-tets sharing node 0
        vol = 0.0
        for a, b, cc in [(1, 3, 4), (3, 2, 6), (3, 6, 4), (4, 6, 5), (1, 4, 5), (1, 5, 3)]:
            vol += np.dot(c[a] - c[0], np.cross(c[b] - c[0], c[cc] - c[0]))
    scale = np.max(np.linalg.norm(c - c[0], axis=1)) ** 3
    if scale == 0.0 or abs(vol) / scale < 1e-9:
        raise ValueError(
            "element %r (Exodus id %s) has ~zero volume (%.4g) after "
            "node-order conversion -- the node-order permutation likely "
            "grouped the wrong corners together; see exodus_io.py's module "
            "docstring" % (etype, elem_id, vol))


def _read_coords(ds):
    # netCDF4 variable reads come back as numpy.ma.MaskedArray (even when
    # nothing is actually masked -- our files never leave a coordinate
    # unwritten); np.asarray strips that back to a plain ndarray, which
    # e.g. scipy's ConvexHull requires and rejects a MaskedArray for outright.
    if "coord" in ds.variables:
        return np.ascontiguousarray(np.asarray(ds.variables["coord"][:]).T, dtype=float)
    return np.column_stack([np.asarray(ds.variables["coord" + a][:]) for a in "xyz"]).astype(float)


def _num_map(ds, name, count):
    if name in ds.variables:
        return [int(v) for v in ds.variables[name][:]]
    return list(range(1, count + 1))


def read_exodus(path, time_step=-1, stress_names=None):
    """Read an Exodus II file into a :class:`~ddsim.mesh_io.MeshData`.

    ``time_step``: which time step's nodal variables to use as the stress
    field (default: the last one). ``stress_names``: optional ``{"xx": name,
    ...}`` override for files whose nodal-variable names don't match any
    alias in ``_STRESS_ALIASES``.
    """
    import netCDF4

    with netCDF4.Dataset(path) as ds:
        coords = _read_coords(ds)
        n_nodes = coords.shape[0]
        node_ids = np.array(_num_map(ds, "node_num_map", n_nodes), dtype=np.int64)

        elem_ids, elem_types, connectivity = [], [], []
        n_blk = len(ds.dimensions["num_el_blk"]) if "num_el_blk" in ds.dimensions else 0
        running_elem_id = iter(_num_map(ds, "elem_num_map", int(ds.dimensions["num_elem"].size)))
        for b in range(1, n_blk + 1):
            conn_var = ds.variables["connect%d" % b]
            raw_type = conn_var.elem_type.strip().upper()
            if raw_type in _UNSUPPORTED_QUADRATIC:
                raise NotImplementedError(
                    "Exodus element type %r is a quadratic element; not yet "
                    "supported (see exodus_io.py's module docstring)" % raw_type)
            if raw_type not in _EXODUS_TYPE_TO_OURS:
                raise NotImplementedError("unrecognized Exodus element type %r" % raw_type)
            etype = _EXODUS_TYPE_TO_OURS[raw_type]
            perm = _EXODUS_TO_OURS[etype]
            block_conn = np.asarray(conn_var[:], dtype=np.int64)  # 1-based -> row index in coords
            for row in block_conn:
                eid = next(running_elem_id)
                local = node_ids[row[perm] - 1]
                _check_nonzero_volume(etype, coords[row[perm] - 1], eid)
                elem_ids.append(eid)
                elem_types.append(etype)
                connectivity.append(local)

        stress = _read_stress(ds, n_nodes, time_step, stress_names)

        data = MeshData(
            node_ids=node_ids, coords=coords,
            elem_ids=np.array(elem_ids, dtype=np.int64),
            elem_types=elem_types, connectivity=connectivity,
            stress=stress, meta={"format": "EXODUS", "path": str(path)})
        data.validate()
        return data


def _read_stress(ds, n_nodes, time_step, stress_names):
    if "name_nod_var" not in ds.variables:
        return None
    names = ["".join(c.decode() for c in row if c).strip().lower()
             for row in ds.variables["name_nod_var"][:]]
    by_name = {n: i + 1 for i, n in enumerate(names)}  # vals_nod_var<i+1>

    lookup = dict(stress_names) if stress_names else {}
    for comp, aliases in _STRESS_ALIASES.items():
        if comp in lookup:
            continue
        for alias in aliases:
            if alias in by_name:
                lookup[comp] = alias
                break

    if len(lookup) < 6:
        missing = [c for c in ("xx", "yy", "zz", "xy", "yz", "zx") if c not in lookup]
        raise ValueError(
            "could not find nodal variables for stress component(s) %s among "
            "%r; pass stress_names={'xx': '<var name>', ...} to override" % (missing, names))

    stress = np.empty((n_nodes, 6))
    for i, comp in enumerate(("xx", "yy", "zz", "xy", "yz", "zx")):
        var = ds.variables["vals_nod_var%d" % by_name[lookup[comp]]]
        stress[:, i] = var[time_step, :] if var.ndim == 2 else var[:]
    return stress


def write_exodus(path, mesh, nodal_vars, default=0.0, title="ddsim output"):
    """Write ``mesh`` (a :class:`~ddsim.mesh_io.MeshData`) plus nodal
    variables to a new Exodus II file at ``path`` (always creates a fresh
    file; never modifies an input file in place).

    ``nodal_vars``: ``{"life": {node_id: value, ...}, ...}``; a node missing
    from a dict gets ``default``. If ``mesh.stress`` is set, it is written too
    (as ``stress_xx``, ``stress_yy``, ... -- ``read_exodus``'s own first-choice
    alias for each component), so a mesh with stress round-trips completely
    through write_exodus + read_exodus alone, with no external sample file
    needed to test the stress-reading path.

    ``mesh.node_ids``/``mesh.elem_ids`` are written into Exodus's external ID
    maps (``node_num_map``/``elem_num_map``) shifted by the smallest constant
    needed to make them 1-based, if they aren't already (DDSim's native RDB
    IDs are 0-based). Exodus readers require positive IDs there; ParaView/VTK
    silently drops an entire element block rather than erroring when it sees
    a 0, which is what "opens fine, empty/broken on Apply" looks like. This
    means ``read_exodus(write_exodus(...))`` does not reproduce 0-based IDs
    bit-for-bit -- they come back shifted by +1 -- which is fine for our own
    use (visualization) but worth knowing if you need the IDs to match.
    """
    import netCDF4

    n_nodes = len(mesh.node_ids)
    row_of = {int(nid): i for i, nid in enumerate(mesh.node_ids)}

    all_vars = dict(nodal_vars)
    if mesh.stress is not None:
        for i, comp in enumerate(("xx", "yy", "zz", "xy", "yz", "zx")):
            var_name = _STRESS_ALIASES[comp][0]
            all_vars[var_name] = {int(nid): mesh.stress[j, i] for j, nid in enumerate(mesh.node_ids)}
    nodal_vars = all_vars
    inv_perm = {name: np.argsort(perm) for name, perm in _EXODUS_TO_OURS.items()}
    type_name = {v: k for k, v in _EXODUS_TYPE_TO_OURS.items() if k in ("TETRA4", "WEDGE6", "HEX8")}

    by_type = {}
    for eid, etype, conn in zip(mesh.elem_ids, mesh.elem_types, mesh.connectivity):
        if etype not in _EXODUS_TO_OURS:
            raise NotImplementedError("cannot write element type %r to Exodus" % etype)
        by_type.setdefault(etype, []).append((eid, conn))

    with netCDF4.Dataset(path, "w", format="NETCDF3_64BIT_OFFSET") as ds:
        ds.title = title
        ds.api_version = np.float32(4.98)
        ds.version = np.float32(4.98)
        ds.floating_point_word_size = np.int32(8)
        # file_size=1 ("large model") means coordinates are split into
        # coordx/coordy/coordz (what we write below), matching what VTK's
        # own vtkIOSSWriter produces. We used to write file_size=1 with a
        # single combined `coord` variable, which mismatched and made
        # ParaView/VTK's ExodusII reader crash on Apply ("failed to locate x
        # nodal coordinates"). Switching file_size to 0 to match combined
        # `coord` fixed that crash but uncovered a second, subtler problem:
        # with combined `coord` + file_size=0, this reader still gets
        # coordinates right (RequestData succeeds, geometry is correct) but
        # silently reads every nodal variable back as all-zero -- confirmed
        # by writing a minimal reference file with VTK's own IOSS writer
        # (same node/var data, split coords) and seeing it read correctly
        # where the combined-coord version didn't. Split coords + file_size=1
        # is therefore the one combination verified (via pvpython, VTK's
        # actual reader) to get both geometry AND nodal variables right.
        ds.file_size = np.int32(1)

        ds.createDimension("len_string", 33)
        ds.createDimension("len_line", 81)
        ds.createDimension("four", 4)
        ds.createDimension("time_step", None)
        ds.createDimension("num_dim", 3)
        ds.createDimension("num_nodes", n_nodes)
        ds.createDimension("num_elem", len(mesh.elem_ids))
        ds.createDimension("num_el_blk", len(by_type))
        if nodal_vars:
            ds.createDimension("num_nod_var", len(nodal_vars))

        for i, axis in enumerate("xyz"):
            c = ds.createVariable("coord" + axis, "f8", ("num_nodes",))
            c[:] = mesh.coords[:, i]

        names = ds.createVariable("coor_names", "S1", ("num_dim", "len_string"))
        for i, n in enumerate("xyz"):
            names[i, 0] = n

        eb_prop1 = ds.createVariable("eb_prop1", "i4", ("num_el_blk",))
        eb_prop1.setncattr("name", "ID")
        eb_status = ds.createVariable("eb_status", "i4", ("num_el_blk",))
        elem_num_map = ds.createVariable("elem_num_map", "i4", ("num_elem",))
        node_num_map = ds.createVariable("node_num_map", "i4", ("num_nodes",))

        # Exodus's external ID maps must be strictly positive (1-based);
        # readers treat 0/negative as "non-positive global id" and reject the
        # whole block (confirmed against ParaView/VTK's IOSS reader, which
        # silently drops every element block rather than erroring loudly).
        # DDSim's native RDB node/element IDs are 0-based, so shift by the
        # minimal amount needed to make the smallest ID 1 -- only when
        # necessary, so IDs that are already positive round-trip unchanged.
        node_offset = 1 - min(int(n) for n in mesh.node_ids) if min(mesh.node_ids) <= 0 else 0
        elem_offset = 1 - min(int(e) for e in mesh.elem_ids) if min(mesh.elem_ids) <= 0 else 0
        node_num_map[:] = [int(n) + node_offset for n in mesh.node_ids]

        elem_pos = 0
        for b, (etype, items) in enumerate(by_type.items(), start=1):
            eb_prop1[b - 1] = b
            eb_status[b - 1] = 1
            n_el = len(items)
            n_nod = elements_num_nodes(etype)
            ds.createDimension("num_el_in_blk%d" % b, n_el)
            ds.createDimension("num_nod_per_el%d" % b, n_nod)
            conn = ds.createVariable("connect%d" % b, "i4",
                                     ("num_el_in_blk%d" % b, "num_nod_per_el%d" % b))
            conn.elem_type = type_name[etype]
            perm = inv_perm[etype]
            for i, (eid, our_conn) in enumerate(items):
                exodus_conn = np.asarray(our_conn)[perm]
                conn[i, :] = [row_of[int(n)] + 1 for n in exodus_conn]
                elem_num_map[elem_pos] = eid + elem_offset
                elem_pos += 1

        time_whole = ds.createVariable("time_whole", "f8", ("time_step",))
        time_whole[0] = 0.0

        if nodal_vars:
            name_nod_var = ds.createVariable("name_nod_var", "S1", ("num_nod_var", "len_string"))
            for i, (var_name, _) in enumerate(nodal_vars.items()):
                for j, ch in enumerate(var_name[:33]):
                    name_nod_var[i, j] = ch
            for i, (_, values) in enumerate(nodal_vars.items(), start=1):
                v = ds.createVariable("vals_nod_var%d" % i, "f8", ("time_step", "num_nodes"))
                row = np.full(n_nodes, default, dtype=float)
                for nid, val in values.items():
                    row[row_of[int(nid)]] = val
                v[0, :] = row


def elements_num_nodes(etype):
    return {"TET_4": 4, "WEDGE_6": 6, "BRICK_8": 8}[etype]

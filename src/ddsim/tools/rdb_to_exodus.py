"""Convert an RDB-format mesh (.con/.nod/.sig/.smp/.edg) to Exodus II, e.g. for
visualization in ParaView -- no crack growth involved, just the mesh and
whatever nodal stress the .sig file carries.

``write_exodus`` writes a MeshData's own ``.stress`` automatically under the
canonical component names (``stress_xx``, ...) whenever it is present, so
this is a thin wrapper: read RDB, write Exodus.

Usage::

    python -m ddsim.tools.rdb_to_exodus <rdb_base> <out.exo>

``<rdb_base>`` is the shared prefix of the .con/.nod/.sig/.smp files (and
.edg, if the mesh has quadratic elements), e.g. ``path/to/SIPS3002`` for
``path/to/SIPS3002.con`` etc.
"""
import sys

from .. import exodus_io, mesh_io


def convert(rdb_base, out_path):
    data = mesh_io.read_rdb(rdb_base)
    exodus_io.write_exodus(out_path, data, nodal_vars={},
                           title="ddsim: %s (converted from RDB)" % rdb_base)
    return data


def main():
    if len(sys.argv) != 3:
        print("usage: python -m ddsim.tools.rdb_to_exodus <rdb_base> <out.exo>")
        sys.exit(1)
    rdb_base, out_path = sys.argv[1], sys.argv[2]
    data = convert(rdb_base, out_path)
    print("wrote %s: %d nodes, %d elements (%s)" % (
        out_path, len(data.node_ids), len(data.elem_ids),
        ", ".join(sorted(set(data.elem_types)))))


if __name__ == "__main__":
    main()

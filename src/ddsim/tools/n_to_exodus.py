"""Convert a real DDSim ``.N`` life-prediction result file (from the actual
2007 parallel Monte Carlo runs, e.g.
``SIPS_data/DDSimLI/SIPS3002_open/ConstantAmplitude/10000_Particles/*.N``)
plus its RDB mesh to Exodus, for visualizing the already-computed results in
ParaView -- no crack growth simulation is run here.

``.N`` file format: one ``doid rid N`` triple per line (whitespace-separated;
``doid`` the node id, ``rid`` a particle/realization id or ``-1`` for a
single deterministic/no-particle-nearby result, ``N`` the predicted life in
cycles). A node can appear on any number of lines, one per Monte Carlo
particle that seeded a flaw there.

Per-node value written to Exodus: the plain arithmetic mean of ``N`` across
every line for that doid -- this is exactly ``Statistic.StatN.SampleMean``
(the reduction ``DamModel.LifeValues`` itself uses for the Monte Carlo case:
plain ``mean(values)``, no weighting), so this reproduces the same per-node
statistic DamMo would have written via ``ToExodusFile``, just from already-
captured output instead of a live run. Nodes absent from the ``.N`` file
(never seeded/run) are written as NaN, matching ``ToExodusFile``'s own
"never run" convention (ParaView shows NaN as masked/blank).

Usage::

    python -m ddsim.tools.n_to_exodus <rdb_base> <n_file> <out.exo> [var_name]

``var_name`` (default ``"life"``) names the nodal variable written to Exodus.
"""
import sys

import numpy as np

from .. import exodus_io, mesh_io


def read_n_file(path):
    """Parse a ``.N`` file into ``{doid: [N, N, ...]}``."""
    per_node = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            doid, _rid, n = int(parts[0]), int(parts[1]), float(parts[2])
            per_node.setdefault(doid, []).append(n)
    return per_node


def node_means(per_node):
    """Per-node arithmetic mean of N, matching Statistic.StatN.SampleMean."""
    return {doid: float(np.mean(values)) for doid, values in per_node.items()}


def convert(rdb_base, n_path, out_path, var_name="life"):
    mesh = mesh_io.read_rdb(rdb_base)
    per_node = read_n_file(n_path)
    values = node_means(per_node)
    exodus_io.write_exodus(
        out_path, mesh, {var_name: values}, default=float("nan"),
        title="ddsim: %s (%s, from %s)" % (rdb_base, var_name, n_path))
    return mesh, values


def main():
    if len(sys.argv) not in (4, 5):
        print("usage: python -m ddsim.tools.n_to_exodus <rdb_base> <n_file> <out.exo> [var_name]")
        sys.exit(1)
    rdb_base, n_path, out_path = sys.argv[1:4]
    var_name = sys.argv[4] if len(sys.argv) == 5 else "life"

    mesh, values = convert(rdb_base, n_path, out_path, var_name)
    vals = list(values.values())
    print("wrote %s: %d/%d nodes have %r (mean %.1f, min %.1f, max %.1f)" % (
        out_path, len(values), len(mesh.node_ids), var_name,
        np.mean(vals), np.min(vals), np.max(vals)))


if __name__ == "__main__":
    main()

"""Pull a nodal variable (default ``"life"``) back out of an Exodus II file
(e.g. the one DDSim now writes by default -- see ``DamMo.ToExodusFile``,
``docs/PORTING_NOTES.md``) as a ``.N``-style text file -- the reverse
direction of ``ddsim.tools.n_to_exodus``, which goes ``.N`` -> Exodus.

Exodus stores one aggregated value per node, with no per-particle
breakdown, so the output always has exactly one line per node, with
``rid=-1`` -- the same "deterministic/no single-particle breakdown"
convention ``NFile``'s own non-Monte-Carlo branch already uses. A node
written as NaN (``ToExodusFile``'s "never run" sentinel) is skipped
entirely, matching a doid that was never in ``doid_list`` at all rather
than one that ran and got an explicit ``-1``.

Usage::

    python -m ddsim.tools.exodus_to_n <exodus_file> <out.N> [var_name]

``var_name`` (default ``"life"``) names the nodal variable to extract.
"""
import sys

from .. import exodus_io


def convert(exodus_path, out_path, var_name="life"):
    values = exodus_io.read_nodal_variable(exodus_path, var_name)
    with open(out_path, "w") as f:
        for doid in sorted(values):
            v = values[doid]
            if v == v:  # skip NaN ("never run" -- see module docstring)
                f.write("%d -1 %r\n" % (doid, v))
    return values


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) not in (2, 3):
        print(__doc__)
        sys.exit(1)
    exodus_path, out_path = argv[0], argv[1]
    var_name = argv[2] if len(argv) == 3 else "life"
    convert(exodus_path, out_path, var_name)


if __name__ == "__main__":
    main()

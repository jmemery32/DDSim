"""Find the N most critical (lowest predicted life) doids from an existing
.N result file, re-run each one's own representative initial flaw through
the CURRENT code, and write each one's crack-growth history as a
-crack_path-style VTK file -- all from a single mesh load, rather than N
separate ``ddsim`` subprocess invocations.

For Monte Carlo (monte=1/2) inputs, "representative" means the doid's own
*largest*-particle initial flaw size -- the one actually simulated cycle-by-
cycle by ``DDSim.Var_Amplitude``/``DamModel.VarAmp`` to determine the doid's
own reported ``Life`` (every *other*, smaller particle at that doid gets an
*interpolated* life estimate from that one curve, never its own simulated
growth history -- see ``DamModel.InterpolateLife`` -- so there is no crack
path to extract for them). This module reproduces exactly that one
simulation per doid directly, instead of running the real driver
(``DDSim.Var_Amplitude``'s monte=1/2 branch *also* re-simulates every
smaller particle in a loop, overwriting ``DamOro[doid].DamEl`` each time via
``Flush()`` -- so by the time it returns, the live history left behind is
whichever particle happened to be simulated *last*, not the largest one;
seen and confirmed while building this tool, see docs/PORTING_NOTES.md).

Usage::

    python -m ddsim.tools.top_crack_paths <rdb_base> <n_file> <par_base> \\
        <conpath> <parpath> <out_dir> [--val <spectrum.val>] [--top N] \\
        [--nore] [--scale S]

Omit ``--val`` for constant-amplitude (RK5) instead of variable-amplitude.
"""
import argparse
import os

from .. import DamMo, MeshTools, Parameters, Statistic, VarAmplitude
from .n_to_exodus import node_means, read_n_file


def rank_critical_doids(n_file_path, top):
    """[(doid, mean_life), ...] for the `top` lowest-life doids, ascending."""
    means = node_means(read_n_file(n_file_path))
    return sorted(means.items(), key=lambda kv: kv[1])[:top]


def _read_particle_filter(parpath, par_base, sets, samples):
    """Same .rnd/.map parsing as DDSim.InitializingStuff's monte==2 branch,
    without that function's sys.argv coupling."""
    ais = Statistic.Stata(sets, samples)
    with open(os.path.join(parpath, par_base + ".rnd")) as f:
        for i in range(sets):
            ais.UpdateSet()
            for _j in range(samples):
                rid, a = f.readline().split()
                ais.UpdateSamples(float(a), i, int(rid))
    ais_map = {}
    with open(os.path.join(parpath, par_base + ".map")) as f:
        for line in f:
            parts = line.split()
            ais_map[int(parts[0])] = [int(x) for x in parts[2:]]
    return ais, ais_map


def _make_local_ais(ais, ais_map, doid):
    """Same logic as DDSim.MakeLocalAis, reproduced to avoid importing
    DDSim.py (which would pull in its module-level argv parsing)."""
    if doid not in ais_map:
        return Statistic.Stata(0, 0)
    local = Statistic.Stata(1, len(ais_map[doid]))
    local.UpdateSet()
    for rid in ais_map[doid]:
        local.UpdateSamples(ais.rid_ai[rid][1], 0, rid)
    return local


def write_top_crack_paths(rdb_base, n_file_path, par_base, conpath, parpath,
                          out_dir, val_path=None, top=50, nore=False,
                          scale=1.0, verbose=True):
    """Returns a list of (doid, historical_mean_life, ai_used, rerun_N,
    willgrow, vtk_path) -- vtk_path is None for a doid with no particles
    mapped (nothing to simulate)."""
    critical = rank_critical_doids(n_file_path, top)
    if verbose:
        print("%d most critical doids (lowest historical mean life):" % len(critical))
        for doid, life in critical:
            print("  doid %-8d life=%.1f" % (doid, life))

    model = MeshTools.MeshTools(conpath + rdb_base, "RDB")
    model.SetPointInsideTolerance(1.0e-7)
    params = Parameters.Parameters(par_base, parpath)

    ais = ais_map = None
    if params.monte in (1, 2):
        ais, ais_map = _read_particle_filter(parpath, par_base, params.sets,
                                             params.samples)

    spec = VarAmplitude.Spectrum(val_path) if val_path else None
    os.makedirs(out_dir, exist_ok=True)

    results = []
    for doid, hist_life in critical:
        # A fresh DamModel per doid -- cheap (no mesh reload) and keeps
        # WriteCrackPathVTK's "pick the first key in DamOro" convention
        # trivially correct without juggling multiple doids in one model.
        cracks = DamMo.DamModel(model, model.GetNodeList(), 0, False, "0",
                                False, params, False, False, "RK5")
        xyz, _delxyz, _sigxyz = model.GetNodeInfo(doid)

        if params.monte in (1, 2):
            local_ais = _make_local_ais(ais, ais_map, doid)
            rv = local_ais.rv
            if not rv or not rv[0]:
                if verbose:
                    print("  doid %d: no particles mapped, skipping" % doid)
                results.append((doid, hist_life, None, None, None, None))
                continue
            ai = max(rv[0].keys())
            cracks.AddDamOro(doid, local_ais)
        else:
            ai = params.a_b[0][0]
            cracks.AddDamOro(doid)

        cracks.AddFDam(doid, xyz, ai, ai, params.material, 0, False)

        if spec is not None:
            N, _fin = cracks.VarAmp(doid, params.Max_crack_size,
                                    params.N_max, spec, nore, scale, params.r)
        else:
            N = cracks.SimDamGrowth(doid, params, "err", scale,
                                    params.material[16])

        vtk_path = os.path.join(out_dir, "doid_%d_crack_path.vtk" % doid)
        cracks.WriteCrackPathVTK(vtk_path, doid)
        willgrow = cracks.DamOro[doid].WillGrow
        if verbose:
            print("  doid %d: historical life=%.1f, rerun (ai=%.4e) "
                  "N=%s willgrow=%s -> %s" %
                  (doid, hist_life, ai, N, willgrow, vtk_path))
        results.append((doid, hist_life, ai, N, willgrow, vtk_path))

    return results


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("rdb_base")
    p.add_argument("n_file")
    p.add_argument("par_base")
    p.add_argument("conpath")
    p.add_argument("parpath")
    p.add_argument("out_dir")
    p.add_argument("--val", default=None, help="variable-amplitude spectrum "
                   "(.val) path; omit for constant-amplitude")
    p.add_argument("--top", type=int, default=50)
    p.add_argument("--nore", action="store_true")
    p.add_argument("--scale", type=float, default=1.0)
    args = p.parse_args(argv)

    conpath = args.conpath if args.conpath.endswith(("/", "\\")) else args.conpath + "/"
    parpath = args.parpath if args.parpath.endswith(("/", "\\")) else args.parpath + "/"

    write_top_crack_paths(args.rdb_base, args.n_file, args.par_base, conpath,
                          parpath, args.out_dir, val_path=args.val,
                          top=args.top, nore=args.nore, scale=args.scale)


if __name__ == "__main__":
    main()

"""In-process multiprocessing replacement for DDSim's old Windows/MPI-cluster
parallel workflow. The original (2007) scheme needed a shared network drive,
`mpirun`, an `MSTI_RANK` env var, and `legacy/windows_cluster_scripts/*.bat`
to partition a doid list across ranks (`node_partition.<rank>` files) and
`ddsim.tools.twins` to merge per-rank output text files back together
afterward. This module does the same job -- partition a doid list, run each
partition's crack-growth simulation, merge the results -- entirely in one
process's memory on a single machine, via `-j <N>` on the `ddsim` CLI.

See docs/PORTING_NOTES.md for why per-doid results can't just be pickled and
returned across the process boundary as-is (DamModel.DamOro[doid] holds a
reference to the whole mesh via its DamEl list), and how
DamMo._DamOroContainer.StripDamElForTransport() works around it.
"""
import random
from concurrent.futures import ProcessPoolExecutor, as_completed

from . import MeshTools, Parameters, DamMo, VarAmplitude

_worker_state = None  # set once per worker process by _worker_init
_worker_state_va = None  # set once per worker process by _worker_init_va


def partition_doids(doid_list, num_workers, seed=None):
    """In-memory equivalent of the old DDSim.DoParallel: randomly shuffle
    doid_list (its only load-balancing heuristic -- "a simple measure for
    speedup", per the original comment) then split into num_workers
    contiguous chunks, each individually sorted (matching the legacy
    node_partition.<rank> files' own convention). Returns a list of up to
    num_workers lists; does not mutate doid_list. Callers should clamp
    num_workers <= len(doid_list) first (unlike the legacy DoParallel, which
    didn't and could produce empty partitions for excess ranks)."""
    if seed:
        random.seed(seed)
    remaining = list(doid_list)
    shuffled = []
    while remaining:
        shuffled.append(remaining.pop(random.randint(0, len(remaining) - 1)))

    chunks = []
    old = 0
    for i in range(num_workers - 1):
        cut = (i + 1) * (len(shuffled) // num_workers)
        chunks.append(sorted(shuffled[old:cut]))
        old = cut
    chunks.append(sorted(shuffled[old:]))
    return chunks


def _worker_init(conpath, filename, parpath, exodus_in, verbose, verify,
                  int_type, ais, ais_map, ncr, scale, errfile_extension,
                  a_b_override=None):
    """ProcessPoolExecutor initializer: runs once per worker process (not
    once per task), so the ~11s real-mesh-load cost is paid once per worker
    and amortized across every doid that worker processes."""
    global _worker_state
    if exodus_in:
        model = MeshTools.MeshTools(exodus_in, 'EXODUS')
    else:
        model = MeshTools.MeshTools(conpath + filename, 'RDB')
    model.SetPointInsideTolerance(1.0e-7)
    parameters = Parameters.Parameters(filename, parpath)
    if a_b_override is not None:
        # Each worker re-reads the .par file independently (see module
        # docstring), so DDSim.main()'s own -ai/-bi override of the parent's
        # Parameters.a_b never reaches here on its own -- re-apply it.
        parameters.a_b = a_b_override
    # DebugGeomUtils/SVIEW always off in workers: those diagnostic dumps
    # fire once at DamModel construction, off the mesh itself, not per doid
    # -- the parent's own unconditionally-built DamModel already handles
    # them; workers doing so too would just be redundant/wasteful.
    cracks = DamMo.DamModel(model, model.GetNodeList(), verbose, verify,
                            '0', False, parameters, False, False, int_type)
    _worker_state = dict(model=model, parameters=parameters, cracks=cracks,
                        ais=ais, ais_map=ais_map, ncr=ncr, scale=scale,
                        int_type=int_type, verbose=verbose, verify=verify,
                        parpath=parpath, filename=filename,
                        errfile_extension=errfile_extension)


def _process_chunk(doid_chunk):
    """Runs in a worker process. Reuses DDSim.MonteSimulation/
    FwdDeterministic directly -- the exact same dispatch
    DDSim.Fwd_Integration's serial loop body already does -- so a -j>=2 run
    is provably running the identical per-doid code as -j 1, just split
    across processes."""
    from . import DDSim as _ddsim  # deferred: avoids a DDSim<->parallel import cycle
    s = _worker_state
    model, cracks, parameters = s['model'], s['cracks'], s['parameters']
    results = {}
    for doid in doid_chunk:
        xyz, _delxyz, _sigxyz = model.GetNodeInfo(doid)
        if parameters.monte in (1, 2):
            _ddsim.MonteSimulation(s['ais'], s['ais_map'], cracks,
                s['parpath'], s['filename'], '0', s['errfile_extension'],
                parameters, doid, xyz, s['verbose'], s['ncr'], s['scale'],
                s['int_type'], s['verify'], parameters.material[16])
        else:
            _ddsim.FwdDeterministic(cracks, doid, xyz, parameters,
                s['verbose'], s['verify'], s['parpath'], s['filename'],
                s['errfile_extension'], s['scale'])
        cracks.DamOro[doid].StripDamElForTransport()
        results[doid] = cracks.DamOro[doid]
        # bound worker memory -- same reasoning as Fwd_Integration's old
        # local_node/data_base branch, just unconditional now.
        del cracks.DamOro[doid]
    return results


def run_parallel(cracks, ais, ais_map, ncr, doid_list, num_workers, conpath,
                 filename, parpath, exodus_in, errfile_extension, scale,
                 int_type, verify, verbose, parameters, saveall, seed=None,
                 a_b_override=None):
    """Parent-side orchestration, called from DDSim.main() in place of
    Fwd_Integration when -j N (N>=2) is given. Mutates cracks.DamOro in
    place with every worker's merged results, then -- if saveall -- writes
    the same .N/.ori/.ai/.af files Fwd_Integration's serial saveall branch
    would have, once in a batch after all workers finish rather than
    streamed per-doid (final file contents are identical either way; the
    one behavior difference is a crash mid-run no longer leaves partial
    output for doids that had already completed)."""
    num_workers = max(1, min(num_workers, len(doid_list)))
    chunks = [c for c in partition_doids(list(doid_list), num_workers, seed) if c]

    with ProcessPoolExecutor(max_workers=num_workers, initializer=_worker_init,
            initargs=(conpath, filename, parpath, exodus_in, verbose, verify,
                      int_type, ais, ais_map, ncr, scale, errfile_extension,
                      a_b_override)) as ex:
        futures = [ex.submit(_process_chunk, chunk) for chunk in chunks]
        for fut in as_completed(futures):
            cracks.DamOro.update(fut.result())

    cracks.RefreshLifeBounds()

    if saveall:
        # One call per doid (not the 'all' batch form, which writes an extra
        # header line, e.g. "Initial_a 0") -- matches Fwd_Integration's
        # serial saveall branch's file format exactly, just written once
        # after all workers finish rather than streamed per-doid.
        doids = sorted(cracks.DamOro)
        nfile = open(parpath + filename + '.N', 'a')
        OrientFile = open(parpath + filename + '.ori', 'a')
        AiFile = open(parpath + filename + '.ai', 'a')
        AfFile = open(parpath + filename + '.af', 'a')
        for doid in doids:
            cracks.NFile(doid, nfile)
            cracks.WriteRotations(OrientFile, doid)
            cracks.WriteInitialAs(AiFile, parameters.N_max, doid)
            cracks.WriteFinalAs(AfFile, doid)
        nfile.close()
        OrientFile.close()
        AiFile.close()
        AfFile.close()


# ---------------------------------------------------------------------------
# -VarAmp -j: same design as above, kept as parallel (pun intended), separate
# functions rather than branching the CA ones, since the two drivers take
# genuinely different arguments (a spectrum + nore, no int_type) and this
# keeps the already-validated CA path completely untouched. Added 2026, see
# docs/PORTING_NOTES.md.
# ---------------------------------------------------------------------------
def _worker_init_va(conpath, filename, parpath, exodus_in, verbose, verify,
                    ais, ais_map, ncr, scale, nore, val_path,
                    a_b_override=None):
    """ProcessPoolExecutor initializer for -VarAmp -j: like _worker_init, but
    also builds this worker's own VarAmplitude.Spectrum (cheap -- just reads
    and parses the .val file) instead of an Int_type-specific DamModel."""
    global _worker_state_va
    if exodus_in:
        model = MeshTools.MeshTools(exodus_in, 'EXODUS')
    else:
        model = MeshTools.MeshTools(conpath + filename, 'RDB')
    model.SetPointInsideTolerance(1.0e-7)
    parameters = Parameters.Parameters(filename, parpath)
    if a_b_override is not None:
        parameters.a_b = a_b_override
    cracks = DamMo.DamModel(model, model.GetNodeList(), verbose, verify,
                            '0', False, parameters, False, False, 'RK5')
    spec = VarAmplitude.Spectrum(val_path)
    _worker_state_va = dict(model=model, parameters=parameters, cracks=cracks,
                           ais=ais, ais_map=ais_map, ncr=ncr, scale=scale,
                           nore=nore, verbose=verbose, verify=verify,
                           spec=spec)


def _process_chunk_va(doid_chunk):
    """Runs in a worker process. Reuses DDSim.VarAmpOneDoid directly -- the
    exact same per-doid body DDSim.Var_Amplitude's serial loop already
    calls -- so a -j>=2 VA run is provably running the identical per-doid
    code as -j 1, just split across processes."""
    from . import DDSim as _ddsim  # deferred: avoids a DDSim<->parallel import cycle
    s = _worker_state_va
    cracks = s['cracks']
    results = {}
    for doid in doid_chunk:
        _ddsim.VarAmpOneDoid(doid, s['model'], cracks, s['parameters'],
            s['spec'], s['verbose'], s['verify'], s['nore'], s['scale'],
            s['ncr'], s['ais'], s['ais_map'])
        cracks.DamOro[doid].StripDamElForTransport()
        results[doid] = cracks.DamOro[doid]
        del cracks.DamOro[doid]  # bound worker memory, same as _process_chunk
    return results


def run_parallel_va(cracks, ais, ais_map, ncr, doid_list, num_workers,
                    conpath, filename, parpath, exodus_in, scale, nore,
                    verify, verbose, parameters, saveall, val_path,
                    seed=None, a_b_override=None):
    """Parent-side orchestration for -VarAmp -j, called from DDSim.main() in
    place of Var_Amplitude when -j N (N>=2) is given. Same merge/-sv-output
    strategy as run_parallel -- see its own docstring."""
    num_workers = max(1, min(num_workers, len(doid_list)))
    chunks = [c for c in partition_doids(list(doid_list), num_workers, seed) if c]

    with ProcessPoolExecutor(max_workers=num_workers, initializer=_worker_init_va,
            initargs=(conpath, filename, parpath, exodus_in, verbose, verify,
                      ais, ais_map, ncr, scale, nore, val_path,
                      a_b_override)) as ex:
        futures = [ex.submit(_process_chunk_va, chunk) for chunk in chunks]
        for fut in as_completed(futures):
            cracks.DamOro.update(fut.result())

    cracks.RefreshLifeBounds()

    if saveall:
        doids = sorted(cracks.DamOro)
        nfile = open(parpath + filename + '.N', 'a')
        OrientFile = open(parpath + filename + '.ori', 'a')
        AiFile = open(parpath + filename + '.ai', 'a')
        AfFile = open(parpath + filename + '.af', 'a')
        for doid in doids:
            cracks.NFile(doid, nfile)
            cracks.WriteRotations(OrientFile, doid)
            cracks.WriteInitialAs(AiFile, parameters.N_max, doid)
            cracks.WriteFinalAs(AfFile, doid)
        nfile.close()
        OrientFile.close()
        AiFile.close()
        AfFile.close()

Usage
=====

See the repository ``README.md`` for installation. This page covers running
the ``ddsim`` command-line driver once installed; ``ddsim -help`` is always
the authoritative, up-to-date list of switches.

A minimal run
--------------

.. code-block:: console

   $ ddsim -base example1 -conpath examples/example1/ -parpath examples/example1/ \
         -v -doid_list 10,0 -scale 100

``-base``/``-conpath``/``-parpath`` locate the mesh (``.con``/``.nod``/``.sig``/
``.smp``/``.edg``) and the ``.par`` material/control-parameter file.
``-doid_list`` is a comma-separated (no spaces) list of node ids to analyze;
``-scale`` multiplies the input stress field (handy for driving a crack to
instability quickly, as the test suite does). ``-v`` prints each node's
predicted life as it's computed.

Deterministic vs. Monte Carlo
-------------------------------

Whether a run is deterministic or Monte Carlo is controlled by the ``.par``
file's ``monte`` keyword (``0`` deterministic, ``1`` sampled from a Weibull
distribution, ``2`` from a microstructural particle-cracking filter -- see
:doc:`theory`), not a command-line flag.

For a deterministic run, ``-ai <a> [-bi <b>]`` overrides the ``.par`` file's
own initial flaw size from the command line, without editing it:

.. code-block:: console

   $ ddsim -base example1 -conpath ./ -parpath ./ -doid_list 10 -scale 100 -ai 0.02

``-bi`` sets the second crack dimension independently; omit it for an
initially circular/equal-dimension flaw. Ignored (with a warning) for
Monte Carlo runs, where ``a_b`` means something different (the sampling
distribution's shape).

Constant vs. variable amplitude
---------------------------------

The default is constant-amplitude loading with adaptive RK5 integration.
``-VarAmp <spectrum file>`` switches to variable-amplitude, cycle-by-cycle
integration instead (see :doc:`theory`); ``-nore`` disables the Willenborg
retardation model for either one.

Running faster: ``-j``
------------------------

.. code-block:: console

   $ ddsim -base SIPS3002 -conpath ... -parpath ... -S -sv -j 4

splits ``doid_list`` across ``4`` worker processes on the local machine
(:mod:`ddsim.parallel`), replacing the original Windows/MPI cluster
workflow -- for constant *or* variable amplitude (``-VarAmp``) alike.
Results are bit-identical to a serial run; see ``docs/PORTING_NOTES.md``
for the design and the real bugs this port found and fixed along the way.
Omit ``-j`` (or pass ``-j 1``) for today's plain serial behavior.

Visualizing results in ParaView
----------------------------------

``-exodus_out <path>`` writes predicted life as an Exodus II nodal variable,
directly viewable as a contour plot in ParaView (``-exodus <path>`` reads
an Exodus mesh/stress field as input, instead of the original ASCII RDB
format). ``ddsim.tools.n_to_exodus`` does the same for an already-computed
``.N`` result file, with no simulation involved.

``-crack_path <path>`` (single doid, no ``-j``) writes that node's full
crack-growth history as a legacy VTK PolyData file -- one polyline per
growth step, colored by cumulative cycle count -- for visualizing the
predicted crack path itself, loaded in ParaView alongside the Exodus mesh:

.. code-block:: console

   $ ddsim -base example1 -conpath ./ -parpath ./ -doid_list 10 -scale 100 \
         -crack_path crack_path.vtk

See :meth:`ddsim.DamMo.DamModel.WriteCrackPathVTK` and
``docs/PORTING_NOTES.md`` for how the per-step cycle count is reconstructed
and why this only works without ``-j``.

To see the *N most critical* crack paths over a whole model at once (e.g.
"the 50 lowest-life nodes") rather than one doid you already picked,
:mod:`ddsim.tools.top_crack_paths` ranks doids by life from an existing
``.N`` result file and writes each one's own representative crack path --
all from a single mesh load:

.. code-block:: console

   $ python -m ddsim.tools.top_crack_paths SIPS3002 historical.N sips3002 \
         ./ ./ ./out_dir --val sips3002.val --top 50

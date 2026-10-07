API Reference
==============

Auto-generated from docstrings. See ``docs/architecture.md`` for a guided
tour of how these modules fit together (the pipeline, and why the numeric
core of the hot-path modules is split into a separate numba-compiled
module); this page is for looking up a specific class or function.

Pipeline
--------

.. autosummary::
   :toctree: generated
   :recursive:

   ddsim.DDSim
   ddsim.Parameters
   ddsim.MeshTools
   ddsim.mesh_io
   ddsim.exodus_io
   ddsim.elements
   ddsim.GeomUtils
   ddsim.DamMo
   ddsim.DamClass
   ddsim.DamHistory
   ddsim.DamErrors
   ddsim.MydadN
   ddsim.Integration
   ddsim.Newton
   ddsim.Statistic
   ddsim.VarAmplitude
   ddsim.parallel

Vector / tensor helpers
------------------------

.. autosummary::
   :toctree: generated
   :recursive:

   ddsim.Vec3D
   ddsim.ColTensor
   ddsim.JohnsVectorTools

Tools
-----

.. autosummary::
   :toctree: generated
   :recursive:

   ddsim.tools.PDF
   ddsim.tools.twins
   ddsim.tools.rdb_to_exodus
   ddsim.tools.n_to_exodus
   ddsim.tools.top_crack_paths

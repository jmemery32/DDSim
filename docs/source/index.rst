DDSim Level I
==============

DDSim Level I is a hierarchical, probabilistic fatigue-life predictor: given
a linear-elastic finite-element stress field, it seeds an elliptical flaw at
each candidate node (or a Monte Carlo set of flaws, from an analytical
distribution or a microstructural particle-cracking filter), grows the
crack under NASGRO/Willenborg fatigue crack growth with constant- or
variable-amplitude loading, and reports the predicted life. It is a
from-scratch Python 3 port of the dissertation codebase described in
Emery et al. (2009) -- see :doc:`theory` for the full citation and the
method's background.

.. toctree::
   :maxdepth: 2

   theory
   usage
   api/index

For the port itself -- what changed from the original 2007 code, what was
verified against it, and every known behavioral difference -- see
``docs/PORTING_NOTES.md`` and ``docs/architecture.md`` in the repository
(not part of this generated site, but linked throughout it).

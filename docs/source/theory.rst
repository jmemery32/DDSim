Theory
======

This page summarizes the fatigue-life prediction method this package
implements and points each piece at the module that implements it. It is a
summary, not a substitute for the source papers in ``docs/papers/`` -- see
the :ref:`references` at the bottom for full citations, and
``docs/PORTING_NOTES.md`` for how (and why) this Python 3 port's behavior
sometimes differs, deliberately, from the originals.

Motivation: damage prognosis, not just damage tolerance
---------------------------------------------------------

Traditional damage-tolerance analysis assumes a flaw already exists (from an
intrinsic material defect, machining, or assembly) and asks how many cycles
it takes to grow to failure, using hand-book stress-intensity-factor (SIF)
solutions and a fast, approximate fatigue code [emery2009]_ --
accurate within the limits of that approximate geometry, but not able to
exploit the detail available in a real finite-element model of a structure,
and typically deterministic rather than reflecting how genuinely random the
underlying material and flaw population actually is.

DDSim was conceived as a step toward *damage prognosis*: predicting a real
structure's remaining useful life by combining physically-grounded,
probabilistic models of damage evolution with an as-built finite-element
model, down to the microstructural scale where damage actually nucleates
[emery2011]_. Level I -- the only level this package implements, and
the only one that was ever matured past an incomplete prototype -- is the
first, fastest, most conservative tier of a planned three-level hierarchy:
an initial, reduced-order screen of an entire structure's finite-element
model to find which locations are actually life-limiting, before spending
any higher-fidelity (Level II/III) analysis effort only where it matters.

The crack-growth model
-----------------------

Crack geometry
^^^^^^^^^^^^^^

Level I idealizes every initial flaw as an elliptical crack, in one of three
analytical forms, chosen from the flaw's location in the finite-element mesh
and its orientation relative to the local principal stresses
[emery2009]_:

* **Fellipse** -- a fully embedded ellipse, using the weight-function
  solution of Roy and Saha for an elliptical crack in an infinite body.
* **Hellipse** -- a semi-elliptical surface crack, using the Raju-Newman
  solution for a surface flaw in a finite-thickness plate under combined
  tension and bending.
* **Qellipse** -- a quarter-elliptical corner crack, using the Raju-Newman
  solution for a corner flaw at a finite-thickness plate's edge.

Which one applies is re-evaluated as the crack grows: an embedded flaw that
breaks through to the surface converts from ``Fellipse`` to ``Hellipse``,
and the angle the local geometry subtends (roughly: quarter-elliptical near
a sharp corner, semi-elliptical as that angle opens back up toward a flat
surface) selects between ``Hellipse``/``Qellipse`` for a surface-originating
flaw. See :class:`ddsim.DamClass.Fellipse`, :class:`ddsim.DamClass.Hellipse`,
:class:`ddsim.DamClass.Qellipse`, and ``DamClass.py``'s own
``GeometryCheck``/``WillGrow`` for the actual, authoritative thresholds --
this summary gives the idea, the code is the ground truth.

The NASGRO growth rate law
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Crack growth rate is computed with the NASGRO equation (Forman-Mettu),
which extends the Paris power law with terms for near-threshold behavior and
the approach to fracture toughness:

.. math::

   \frac{da}{dN} = C (\Delta K_{\mathrm{eff}})^n
     \frac{\left(1 - \dfrac{\Delta K_{th}}{\Delta K}\right)^p}
          {\left(1 - \dfrac{K_{max}}{K_c}\right)^q}

where :math:`\Delta K_{\mathrm{eff}} = \Delta K \left(\dfrac{1-f}{1-R}\right)`
folds in Newman's crack-opening (closure) model, :math:`R = K_{min}/K_{max}`,
:math:`K_c` is the critical fracture toughness, and :math:`C`, :math:`n`,
:math:`p`, :math:`q` are material constants [emery2009]_. This is
implemented in :class:`ddsim.MydadN.dadN`.

Level I is applied to flaws far smaller than the large-crack data this
equation was fit to, where growth is dominated by microstructural
randomness rather than plasticity-induced closure, and tends to grow
*faster* than Eq. above predicts near threshold. Rather than use linear
elastic fracture mechanics at a length scale where it's genuinely dubious,
Level I instead defines an *effective load ratio* that decays from a
conservative, closure-free value back to the real applied :math:`R` as the
crack accumulates cycles:

.. math::

   R_{\mathrm{eff}} = (1 - e^{-N/\mathtt{damp}})(R - 0.7) + 0.7

so that a young crack (small :math:`N`) sees a growth rate close to
:math:`\Delta K = K_{max} - K_{min}` with no closure benefit, converging
smoothly to the standard NASGRO equation as it ages [emery2009]_. The
``damp`` material parameter in every ``.par`` file controls how quickly that
convergence happens -- see :class:`ddsim.MydadN.dadN` and
``docs/PORTING_NOTES.md`` (the "SIPS3002 model sanity check" section) for
how this exact parameter was investigated and confirmed against the real
SIPS3002 dataset's own ``.par`` files during this port.

Retardation (load interaction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An overload transiently slows subsequent growth by leaving a larger plastic
zone ahead of the crack tip than the current load alone would produce.
Level I accounts for this with the Generalized Willenborg model
[emery2009]_: Willenborg's plastic-zone-based effective stress range
is reformulated directly in terms of the stress-intensity factor, which
lowers the effective :math:`R` for the same nominal :math:`\Delta K` until
the overload's effect decays away. Implemented in
:class:`ddsim.MydadN.Willenborg`; disabled with ``-nore`` on the command
line (constant-amplitude NASGRO growth with no retardation effects).

Integration
^^^^^^^^^^^

For constant-amplitude loading, the cycles-to-grow-by-:math:`da` relation
above is integrated forward in :math:`N` with an adaptive embedded
Runge-Kutta 5(4) (Cash-Karp) scheme, falling back to forward Euler or
classical RK4 on request (:mod:`ddsim.Integration`). This machinery is this
project's own numerical implementation, not something the source papers
specify in this much detail -- see ``docs/PORTING_NOTES.md``'s "RK5
near-instability overshoot" section for real bugs found and fixed in it
during this port, and for why an adaptive RK scheme needs real care right
at the stiff region approaching instability. Variable-amplitude (spectrum)
loading instead integrates cycle-by-cycle (:class:`ddsim.DamClass.Damage`'s
``VarAmp``) -- deliberately *not* through the same RK machinery, since an
arbitrary load spectrum isn't smooth/differentiable in the way an adaptive
continuous integrator assumes.

Stochastic initial flaw size
-----------------------------

Level I can seed each candidate node with a flaw size in three ways
(:class:`ddsim.Parameters.Parameters`, :func:`ddsim.DDSim.Monte`):

* **Deterministic** -- a single, user-specified initial size (the ``.par``
  file's ``a_b``, or the command line's ``-ai``/``-bi``).
* **A three-parameter Weibull distribution** sampled directly,

  .. math::

     f_X(x) = \frac{\rho}{x_0}\left(\frac{x-\gamma}{x_0}\right)^{\rho-1}
              e^{-\left(\frac{x-\gamma}{x_0}\right)^{\rho}}

  with shape :math:`\rho`, scale :math:`x_0`, and truncation :math:`\gamma`
  [emery2009]_, intended to reflect an observed distribution of
  metallurgical flaws (such as second-phase particle sizes) when no more
  specific microstructural model is available.
* **A microstructural particle-cracking filter** -- rather than guessing a
  generic flaw-size distribution, a companion modeling effort
  ([bozek2008]_) predicts, from a statistically representative
  population of constituent particles (size, aspect ratio, and the
  crystallographic orientation of the surrounding grain), *which* particles
  actually crack under a given local strain level, and the resulting
  predicted crack sizes become the initial-flaw population Level I samples
  from (the ``.rnd``/``.map`` "particle filter" files read by
  :func:`ddsim.DDSim.MonteSimulation`). This is the physically-motivated
  alternative to the generic Weibull fit above, for the one alloy
  (7075-T651) it was developed for.

.. _references:

References
----------

The full PDFs are in ``docs/papers/``.

.. [emery2009] J.M. Emery, J.D. Hochhalter, P.A. Wawrzynek, G. Heber,
   A.R. Ingraffea, "DDSim: A hierarchical, probabilistic, multiscale damage
   and durability simulation system -- Part I: Methodology and Level I,"
   *Engineering Fracture Mechanics* 76 (2009) 1500-1530.
   ``docs/papers/Emery, Hochhalter, Wawrzynek & Ingraffea - DDSim Level I.pdf``.
   This is the paper this package implements -- the primary reference for
   everything on this page.

.. [emery2011] J.M. Emery, A.R. Ingraffea, "DDSim: Framework for Multiscale
   Structural Prognosis," in S. Ghosh, D. Dimiduk (eds.), *Computational
   Methods for Microstructure-Property Relationships*, Springer, 2011.
   ``docs/papers/Emery_Ingraffea_DDSim_CompMethodsStructPropRelations.pdf``
   (the published chapter) and ``docs/papers/emery_ingraffea13.pdf`` (the
   authors' own manuscript copy of the same chapter).

.. [bozek2008] J.E. Bozek, J.D. Hochhalter, M.G. Veilleux, M. Liu, G. Heber,
   S.D. Sintay, A.D. Rollett, D.J. Littlewood, A.M. Maniatty, H. Weiland,
   R.J. Christ Jr., J. Payne, G. Welsh, D.G. Harlow, P.A. Wawrzynek,
   A.R. Ingraffea, "A geometric approach to modeling microstructurally small
   fatigue crack formation, part I: probabilistic simulation of constituent
   particle cracking in AA 7075-T651."
   ``docs/papers/Probabilistic simulation of constituent particle cracking.pdf``.

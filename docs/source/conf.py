"""Sphinx configuration for the DDSim documentation.

Build with ``make docs`` (from the repo root) or directly:

    sphinx-build -b html docs/source docs/build/html

then open docs/build/html/index.html.
"""
import os
import sys

# So autodoc/autosummary can import ddsim without installing it first.
sys.path.insert(0, os.path.abspath("../../src"))

project = "DDSim"
copyright = "2026, John M. Emery"
author = "John M. Emery"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
]

# One stub page per module, auto-generated from the autosummary directives
# in api/index.rst -- no hand-maintained per-module .rst files to go stale.
autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
# The source is a straight port of 2007 Python 2 code with plain, informal
# docstrings (not Google/NumPy style) -- napoleon is harmless here but add
# it anyway for any future docstrings that do use those sections.
napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_title = "DDSim Level I"

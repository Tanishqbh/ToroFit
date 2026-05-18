# Configuration file for the Sphinx documentation builder.
#
# For a full list of options see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------

project = "ToroFit"
copyright = "2024, ToroFit Contributors"
author = "ToroFit Contributors"
release = "0.1.0"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "nbsphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

# Re-execute a notebook at build time only if it ships with no outputs.
# Notebooks that can't be reproduced in CI should be committed with their
# outputs preserved (see CONTRIBUTING / installation docs for the
# `keep_output` metadata workflow with nbstripout).
nbsphinx_execute = "auto"

# -- Options for HTML output -------------------------------------------------

html_theme = "alabaster"
html_static_path = ["_static"]

""" Sphinx configuration """

import os
import sys

# Code
sys.path.insert(0, os.path.abspath(".."))

# Project information
project = "super_material"
author = "Paul le Roux <pleroux0@outlook.com>"
release = "1.0.0"

# Extentions
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.imgmath",
    "sphinx_rtd_theme",
    "sphinxcontrib.bibtex",
]

# Configuration
add_module_names = False

# Ignore
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Theme
html_theme = "sphinx_rtd_theme"

# Static content
html_static_path = ["static"]

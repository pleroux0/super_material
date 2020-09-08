""" Sphinx configuration """

import os
import sys

# Code
sys.path.insert(0, os.path.abspath(".."))

# Project information
project = "super_material"
author = "Paul le Roux <pleroux0@outlook.com>"
release = "1.0.0a"

# Extentions
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.githubpages",
    "sphinx.ext.imgmath",
    "sphinx.ext.inheritance_diagram",
    "sphinx_rtd_theme",
    "sphinxcontrib.bibtex",
]

# Configuration
add_module_names = False
autodoc_typehints = "description"
autoclass_content = "both"

# Ignore
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Theme
html_theme = "sphinx_rtd_theme"

# Static content
html_static_path = ["static"]

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import re


# -- Project information -----------------------------------------------------

project = "MIPTools"
copyright = "2021, Bailey Lab"
author = "Bailey Lab"
version = "v0.4.0.9000"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.githubpages",  # add .nojekyll to gh-pages
    "myst_parser",  # write docs using MyST (a flavor of markdown)
    "sphinx_copybutton",  # add copy button to code chunks
    "sphinx_toolbox.github",  # link to github
    "sphinx_licenseinfo",  # add license information
    "notfound.extension",  # 404 page
    "sphinx.ext.autosectionlabel",  # reference sections using their title
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "README.md"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_context = {
    "display_github": True,
    "github_user": "bailey-lab",
    "github_repo": "MIPTools",
    "github_version": "master/docs/",
}


# -- Sphinx Toolbox configuration-----------------------------------------------
github_username = "bailey-lab"
github_repository = "MIPTools"


# -- 404 Page configuration-----------------------------------------------------
notfound_urls_prefix = "/MIPTools/"


# -- Auto Section configuration-------------------------------------------------
# Make sure the target is unique
autosectionlabel_prefix_document = True


# -- Variables available in all rst files---------------------------------------

stable_version = re.search(r"(v\d+.\d+.\d+)", version).group(1)
if stable_version == version:
    container = f"miptools_{stable_version}.sif"
else:
    container = "miptools_dev.sif"

rst_prolog = """
.. |stable_version| replace:: {0}
.. |container| replace:: {1}
""".format(
    stable_version, container
)

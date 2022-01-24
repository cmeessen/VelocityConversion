# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from VelocityConversion import __version__
sys.path.insert(0, os.path.abspath('../../..'))


# -- Project information -----------------------------------------------------

project = 'VelocityConversion'
copyright = '2019, Christian Meeßen'
author = 'Christian Meeßen'

# The full version, including alpha/beta/rc tags
version = __version__
release = version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.githubpages',
    'sphinx.ext.doctest',
    'IPython.sphinxext.ipython_console_highlighting',
    'IPython.sphinxext.ipython_directive',
    'matplotlib.sphinxext.plot_directive',
    'm2r2',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

# -- sphinx_rtd_theme --------------------------------------------------------

html_theme_options = {
    'canonical_url': '',
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    # 'vcs_pageview_mode': '',
    # 'style_nav_header_background': 'white',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

html_context = {
    'display_github': True,
    'github_user': 'cmeessen',
    'github_repo': 'pyGMS',
    'github_version': 'master/docs/sphinx',
    'menu_links_name': 'More',
    'menu_links': [
        (
            '<i class="fa fa-github fw"></i> GitHub repository',
            'https://github.com/cmeessen/VelocityConversion'
        ),
        (
            '<i class="fa fa-bug fw"></i> Open an issue',
            'https://github.com/cmeessen/VelocityConversion/issues/new'
        ),
        (
            '<i class="fa fa-quote-left fw"></i> Cite VelocityConversion',
            'https://www.zenodo.org/badge/latestdoi/87794116'
        )
    ]
}

# -- sphinx.ext.napoleon -----------------------------------------------------
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_admonition_for_examples = False
napoleon_use_ivar = True


# -- copybutton --------------------------------------------------------------
def setup(app):
    """Adds copy buttons to code blocks, code location in ./_static"""
    app.add_css_file('copy_button.css')
    app.add_js_file('copy_button.js')
    app.add_js_file('clipboard.js')

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
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# on_rtd is whether we are on readthedocs.org, this line of code grabbed from docs.readthedocs.org
import os
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]



# -- Project information -----------------------------------------------------
project = 'EnDAS'
copyright = '2020, Martin Gunia'
author = 'Martin Gunia'

# The short X.Y version.
version = '0.2'

# The full version, including alpha/beta/rc tags
release = '0.2'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [ 'sphinx.ext.autodoc',
               'sphinx.ext.autosummary',
               'sphinx.ext.mathjax',
               #"sphinx.ext.imgmath",
               'sphinx.ext.napoleon',
               'breathe',
               'exhale'
]



breathe_projects = {
    "EnDAS": "./doxyoutput/xml"
}

breathe_default_project = "EnDAS"

doxygen_input = '''
INPUT = ../include
JAVADOC_AUTOBRIEF = YES
STRIP_FROM_INC_PATH = ../include
EXTRACT_LOCAL_CLASSES = NO
HIDE_UNDOC_MEMBERS = YES
HIDE_UNDOC_CLASSES = YES
HIDE_FRIEND_COMPOUNDS = YES
SHOW_USED_FILES = NO
'''

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./cppapi",
    "rootFileName":          "cppapi_root.rst",
    "rootFileTitle":         "C/C++ API",
    "doxygenStripFromPath":  "..",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": doxygen_input
}


# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'





autodoc_member_order = "bysource"
autosummary_generate = True

#autoclass_content = "both"

autodoc_default_options = {
    'member-order': 'bysource',
    #'special-members': '__init__',
    #'undoc-members': False,
    'show-inheritance': True,
    #'inherited-members' : True
}

autodoc_inherit_docstrings = False



# Add any paths that contain templates here, relative to this directory.
templates_path = ['templates']

# Specify index.rst as the master file explicitly as it changed in Sphinx 2.0
master_doc = 'index'

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.
# Please note that since version 1.4.0 the sphinx_rtd_theme is no longer
# bundled with Sphinx. Install it with `pip install sphinx_rtd_theme`.
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['static']

html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'top',
    'style_external_links': False,
    # Toc options
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 3,
    'includehidden': True,
    'titles_only': False
}


def setup(app):
   app.add_stylesheet('css/custom.css')

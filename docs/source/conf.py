# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath("../../"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'HIPPO'
copyright = '2024, Max Winokan'
author = 'Max Winokan'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
              'sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.extlinks',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon',
              'sphinx.ext.intersphinx',
              'sphinxcontrib.prettyspecialmethods',
              ]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_logo = '../../logos/hippo_logo_tightcrop.png'
html_favicon = '../../logos/hippo_assets-02.gif'

html_theme_options = {
    'navigation_depth': -1,
    'logo_only': True,
    'prev_next_buttons_location': 'both',
}

html_css_files = [
    'css/custom.css'
]

source_suffix = '.rst'
master_doc = 'index'
# project = 'ASE'
# copyright = f'{datetime.date.today().year}, ASE-developers'
# templates_path = ['templates']
exclude_patterns = ['build']
# default_role = 'math'
# pygments_style = 'sphinx'
autoclass_content = 'both'
modindex_common_prefix = ['hippo.']

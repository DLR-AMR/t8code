# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import subprocess, os

project = 't8code'
copyright = '2024, Johannes Holke, David Knapp, Sandro Elsweijer, Ioannis Lilikakis, Lukas Dreyer, Jakob Fußbroich, Carsten Burstedde, Chiara Hergl, Johannes Markert, Niklas Boeing, Florian Becker, Prasanna Ponnusamy'
author = 'Johannes Holke, David Knapp, Sandro Elsweijer, Ioannis Lilikakis, Lukas Dreyer, Jakob Fußbroich, Carsten Burstedde, Chiara Hergl, Johannes Markert, Niklas Boeing, Florian Becker, Prasanna Ponnusamy'
version = '4.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [ "breathe", "exhale", "sphinx.ext.mathjax", "sphinx.ext.graphviz" ]

templates_path = ['_templates']
exclude_patterns = []


def configureDoxyfile(input_dir, output_dir):
    with open('Doxyfile.in', 'r') as file :
        filedata = file.read()

    filedata = filedata.replace('@DOXYGEN_INPUT_DIR@', input_dir)
    filedata = filedata.replace('@DOXYGEN_OUTPUT_DIR@', output_dir)

    with open('Doxyfile', 'w') as file:
        file.write(filedata)

# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

breathe_projects = {}

if read_the_docs_build:
    input_dir = '../../src'
    output_dir = 'build'
    configureDoxyfile(input_dir, output_dir)
    subprocess.call('doxygen', shell=True)
    breathe_projects["T8code"] = output_dir + '/xml'


html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '../../t8code_logo.png'

# Breathe Configuration
breathe_default_project = "T8code"
breathe_implementation_filename_extensions = ['.c', '.cpp']

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "doxygenStripFromPath":  "/localdata1/knap_da/t8code/t8code",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle":         "Library API",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    "INPUT = ../../src"
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

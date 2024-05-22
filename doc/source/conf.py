# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 't8code'
copyright = '2024, Johannes Holke, David Knapp, Sandro Elsweijer, Ioannis Lilikakis, Lukas Dreyer, Jakob Fußbroich, Carsten Burstedde, Chiara Hergl, Johannes Markert, Niklas Boeing, Florian Becker, Prasanna Ponnusamy'
author = 'Johannes Holke, David Knapp, Sandro Elsweijer, Ioannis Lilikakis, Lukas Dreyer, Jakob Fußbroich, Carsten Burstedde, Chiara Hergl, Johannes Markert, Niklas Boeing, Florian Becker, Prasanna Ponnusamy'
release = '2.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [ "breathe" ]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Breathe Configuration
#sys.path.append("/localdata1/knap_da/projects/t8code/t8code/")
breathe_projects = {"T8code" : "/localdata1/knap_da/projects/t8code/t8code_cmake/doc/xml"}
breathe_default_project = "T8code"
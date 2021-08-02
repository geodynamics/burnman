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
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('../examples'))

import burnman
import burnman.version

# -- Project information -----------------------------------------------------

project = 'BurnMan'
copyright = '2021, Sanne Cottaar, Timo Heister, Robert Myhill, Ian Rose, Cayman Unterborn'
author = 'Sanne Cottaar, Timo Heister, Robert Myhill, Ian Rose, Cayman Unterborn'

# The short X.Y version.
version = burnman.version.short_version
# The full version, including alpha/beta/rc tags
release = burnman.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary',
              'sphinx.ext.mathjax', 'numpydoc', 'sphinx.ext.intersphinx',
              'sphinx.ext.ifconfig', 'sphinx.ext.inheritance_diagram',
              'sphinx.ext.viewcode', 'sphinxcontrib.bibtex']
bibtex_bibfiles = ['ref.bib']

numpydoc_show_class_members = False

autodoc_default_options = {'members': True,
                           'undoc-members': True,
                           'show-inheritance': True,
                           'inherited-members': True}

autodoc_member_order = 'bysource'

# We can mock import modules if they will break the build
# For example, the pycddlib library depends on C modules
# which can't be installed by readthedocs.
autodoc_mock_imports = ['pycddlib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Output file base name for HTML help builder.
htmlhelp_basename = 'BurnMandoc'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Options for LaTeX output --------------------------------------------

preamble1 = '''
\\newcommand{{\\burnmanversion}}{{{0}}}
'''.format(release)
preamble = preamble1 + r'''
\usepackage{textpos}
\usepackage{amssymb}

\newcommand{\burnman}{\texttt{\bf BurnMan}}

  '''

latex_maketitle = r'''
\begin{titlepage}

{

\definecolor{dark_grey}{gray}{0.3}
\definecolor{aspect_blue}{rgb}{0.9,0.35,0.3}

%LINE 1%
{
\renewcommand{\familydefault}{\sfdefault}

\pagenumbering{gobble}
\begin{center}
\resizebox{\textwidth}{!}{\textcolor{dark_grey}{\fontfamily{\sfdefault}\selectfont
COMPUTATIONAL INFRASTRUCTURE FOR GEODYNAMICS (CIG)
}}

\hrule

%LINE 2%
\color{dark_grey}
\rule{\textwidth}{2pt}

%LINE 3%
\color{dark_grey}
% FILL: additional organizations
% e.g.: {\Large Organization 1\\Organization 2}
{\Large }
\end{center}

%COLOR AND CODENAME BLOCK%
\begin{center}
\resizebox{\textwidth}{!}{\colorbox
% FILL: color of code name text box
% e.g. blue
{aspect_blue}{\fontfamily{\rmdefault}\selectfont \textcolor{white} {
% FILL: name of the code
% You may want to add \hspace to both sides of the codename to better center it, such as:
% \newcommand{\codename}{\hspace{0.1in}CodeName\hspace{0.1in}}
\hspace{0.1in}\burnman{}\hspace{0.1in}
}}}
\\[12pt]
{\Large a thermodynamics and thermoelasticity toolkit}
\end{center}

%MAIN PICTURE%
\begin{textblock*}{5.in}(0.3in,0.3in)
% FILL: image height
% e.g. height=6.5in
\begin{center}
\vspace{.5in}
\includegraphics[height=4.5in]{../../burnjack-small.png}
% FILL: image file name
% e.g. cover_image.png
\end{center}
\end{textblock*}

%USER MANUAL%
\color{dark_grey}
\hfill{\Huge \fontfamily{\sfdefault}\selectfont User Manual \\
% FILL: manual version
% e.g. 1.0
\raggedleft \huge \fontfamily{\sfdefault}\selectfont Version {\burnmanversion}\\
%\\\large(generated from subversion: $Revision: 2568 $)\\
}

%AUTHOR(S) & WEBSITE%
\null
\vfill
\color{dark_grey}
\Large \hfill {\raggedleft \fontfamily{\sfdefault}\selectfont
% FILL: author list
% e.g. Author One\\Author Two\\Author Three\\
% be sure to have a newline (\\) after the final author
Sanne Cottaar\\Timo Heister\\Robert Myhill\\Ian Rose\\Cayman Unterborn\\
}

{\fontfamily{\sfdefault}\selectfont \href{http://geodynamics.org}{http://geodynamics.org}}


%LINE%
\color{dark_grey}
\rule{\textwidth}{2pt}

}

\pagebreak
}
\end{titlepage}
'''

language='en'

latex_logo='burnjack-small.png'

latex_elements = {
    'sphinxsetup': '',
    'passoptionstopackages': r'\PassOptionsToPackage{table}{xcolor}',
    # The paper size ('letterpaper' or 'a4paper').
    'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    'pointsize': '11pt',
    'preamble': preamble,
    'maketitle': latex_maketitle
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
    ('index_pdf', 'BurnMan.tex', u'BurnMan Documentation', author, 'manual'),
]


# The name of an image file (relative to this directory) to place at the top of
# the title page.
# latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
# latex_use_parts = False

# If true, show page references after internal links.
# latex_show_pagerefs = False

# If true, show URL addresses after external links.
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
# latex_appendices = []

# If false, no module index is generated.
latex_domain_indices = False


# -- Options for manual page output --------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'burnman', u'BurnMan Documentation', [author], 1)
]

# If true, show URL addresses after external links.
# man_show_urls = False


# -- Options for Texinfo output ------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    ('index', 'BurnMan', u'BurnMan Documentation',
     author, 'BurnMan', 'One line description of project.',
     'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
# texinfo_appendices = []

# If false, no module index is generated.
# texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
# texinfo_show_urls = 'footnote'

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'http://docs.python.org/': None}

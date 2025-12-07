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
import pkgutil
import inspect

# Define the canonical URL if you are using a custom domain on Read the Docs
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")

# Tell Jinja2 templates the build is running on Read the Docs
if os.environ.get("READTHEDOCS", "") == "True":
    if "html_context" not in globals():
        html_context = {}
    html_context["READTHEDOCS"] = True


sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../examples"))
sys.path.insert(0, os.path.abspath("../tutorial"))

import burnman
from burnman.eos.helper import eos_methods, eos_names
from burnman.eos.property_modifiers import modifier_functions, modifier_names

import datetime

today = datetime.date.today()
year = today.year

today = repr(today)
year = str(year)

# -- Project information -----------------------------------------------------

project = "BurnMan"
copyright = (
    f"{year}, Robert Myhill, Sanne Cottaar, Timo Heister, Ian Rose, Cayman Unterborn"
)
author = "Robert Myhill, Sanne Cottaar, Timo Heister, Ian Rose, Cayman Unterborn"

# The short X.Y version.
v = burnman.__version__.split(".")
version = f"{v[0]}.{v[1]}"
# The full version, including alpha/beta/rc tags
release = burnman.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "numpydoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.ifconfig",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.viewcode",
    "sphinxcontrib.bibtex",
    "nbsphinx",
]

autosummary_generate = True

templates_path = ['_templates']

bibtex_bibfiles = ["ref.bib"]

numpydoc_show_class_members = False

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "inherited-members": True,
}

autodoc_member_order = "bysource"

# We can mock import modules if they will break the build
# For example, the pycddlib library depends on C modules
# which can't be installed by readthedocs.
autodoc_mock_imports = ["pycddlib"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix of source filenames.
source_suffix = ".rst"

# The encoding of source files.
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = "index"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

# Output file base name for HTML help builder.
htmlhelp_basename = "BurnMandoc"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# -- Options for LaTeX output --------------------------------------------

preamble1 = """
\\newcommand{{\\burnmanversion}}{{{0}}}
""".format(
    release
)
preamble = (
    preamble1
    + r"""
\usepackage{textpos}
\usepackage{amssymb}

\newcommand{\burnman}{\texttt{\bf BurnMan}}

  """
)

latex_maketitle = r"""
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
\includegraphics[height=4.5in]{burnjack-small.png}
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
Robert Myhill\\Sanne Cottaar\\Timo Heister\\Ian Rose\\Cayman Unterborn\\
}

{\fontfamily{\sfdefault}\selectfont \href{http://geodynamics.org}{http://geodynamics.org}}


%LINE%
\color{dark_grey}
\rule{\textwidth}{2pt}

}

\pagebreak
}
\end{titlepage}
"""

language = "en"

latex_logo = "burnjack-small.png"
html_logo = "burnjack-small.png"

latex_elements = {
    "sphinxsetup": "",
    "passoptionstopackages": r"\PassOptionsToPackage{table}{xcolor}",
    # The paper size ('letterpaper' or 'a4paper').
    "papersize": "letterpaper",
    # The font size ('10pt', '11pt' or '12pt').
    "pointsize": "11pt",
    "preamble": preamble,
    "maketitle": latex_maketitle,
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
    ("index_pdf", "BurnMan.tex", "BurnMan Documentation", author, "manual"),
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
man_pages = [("index", "burnman", "BurnMan Documentation", [author], 1)]

# If true, show URL addresses after external links.
# man_show_urls = False


# -- Options for Texinfo output ------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        "index",
        "BurnMan",
        "BurnMan Documentation",
        author,
        "BurnMan",
        "One line description of project.",
        "Miscellaneous",
    ),
]

# Documents to append as an appendix to all manuals.
# texinfo_appendices = []

# If false, no module index is generated.
# texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
# texinfo_show_urls = 'footnote'

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'python': ('https://docs.python.org/3', None)}

inheritance_graph_attrs = dict(rankdir="TB", size='""')


# Auto-generate calibrant documentation
output_rst = os.path.join(os.path.dirname(__file__), "autogenerated/calibrants.rst")

with open(output_rst, "w") as f:
    # Walk through all submodules of burnman.calibrants
    for loader, name, ispkg in pkgutil.iter_modules(burnman.calibrants.__path__, burnman.calibrants.__name__ + '.'):
        f.write(f".. automodule:: {name}\n")
        f.write("   :no-index:\n")
        f.write("   :no-inherited-members:\n\n")

# Auto-generate mineral documentation
output_rst = os.path.join(os.path.dirname(__file__), "autogenerated/minerals.rst")

with open(output_rst, "w") as f:
    # Walk through all submodules of burnman.minerals
    for loader, name, ispkg in pkgutil.iter_modules(burnman.minerals.__path__, burnman.minerals.__name__ + '.'):
        f.write(f".. automodule:: {name}\n")
        f.write("   :no-index:\n")
        f.write("   :no-inherited-members:\n\n")

# Auto-generate equation of state documentation
output_rst = os.path.join(os.path.dirname(__file__), "autogenerated/equations_of_state.rst")

thermal_equations_of_state = []
isothermal_equations_of_state = []
other_classes = []
functions = []

for loader, name, ispkg in pkgutil.iter_modules(burnman.eos.__path__, burnman.eos.__name__ + '.'):

    if name == "burnman.eos.equation_of_state":
        continue
    if name == "burnman.eos.helper":
        continue

    # For each module, import it and get all the classes defined in it
    module = __import__(name, fromlist=["dummy"])

    for attr in dir(module):
        obj = getattr(module, attr)

        # Only consider objects defined in this module
        if hasattr(obj, "__module__") and obj.__module__ == name:
            # Only consider classes
            if isinstance(obj, type):

                # Check the if the class is derived from burnman.eos.equation_of_state.IsothermalEquationOfState
                if issubclass(obj, burnman.eos.equation_of_state.IsothermalEquationOfState):
                    isothermal_equations_of_state.append((name, attr))

                elif issubclass(obj, burnman.eos.equation_of_state.EquationOfState):
                    thermal_equations_of_state.append((name, attr))

                else:
                    other_classes.append((name, attr))

            # Check if it's a non-hidden function
            elif callable(obj) and attr[0] != "_":
                functions.append((name, attr))

# convert to dicts, where keys are module names and values are lists of class/function names
isothermal_equations_of_state_dict = {}
for mod, cls in isothermal_equations_of_state:
    if mod not in isothermal_equations_of_state_dict:
        isothermal_equations_of_state_dict[mod] = []
    isothermal_equations_of_state_dict[mod].append(cls)

thermal_equations_of_state_dict = {}
for mod, cls in thermal_equations_of_state:   
    if mod not in thermal_equations_of_state_dict:
        thermal_equations_of_state_dict[mod] = []
    thermal_equations_of_state_dict[mod].append(cls)

other_classes_dict = {}
for mod, cls in other_classes:   
    if mod not in other_classes_dict:
        other_classes_dict[mod] = []
    other_classes_dict[mod].append(cls)

functions_dict = {}
for mod, func in functions:
    if mod not in functions_dict:
        functions_dict[mod] = []
    functions_dict[mod].append(func)

output_rst = os.path.join(
    os.path.dirname(__file__), "autogenerated/api_eos.rst"
)

with open(output_rst, "w") as f:

    f.write(".. _ref-api-eos:\n\n")
    f.write("Equation of State Classes and Functions\n")
    f.write("=======================================\n\n")
    f.write(".. contents::\n")
    f.write("   :local:\n")
    f.write("   :class: this-will-duplicate-information-and-it-is-still-useful-here\n")
    f.write("   :depth: 2\n\n\n")
    for i, (title, d) in enumerate([
        ("Isothermal Equations of State", isothermal_equations_of_state_dict),
        ("Thermal Equations of State", thermal_equations_of_state_dict),
        ("Other Equation of State classes", other_classes_dict),
        ("Equation of State Functions", functions_dict),
    ]):
        f.write(f"{title}\n{'-' * len(title)}\n\n")

        if i == 0:
            f.write("Base class\n")
            f.write("~~~~~~~~~~\n\n")
            f.write(".. autoclass:: burnman.eos.IsothermalEquationOfState\n\n")
        elif i == 1:
            f.write("Base class\n")
            f.write("~~~~~~~~~~\n\n")
            f.write(".. autoclass:: burnman.eos.EquationOfState\n\n")

        for mod, items in d.items():
            name = mod.split(".")[-1].replace("_", " ")
            f.write(f"{name}\n{'~' * len(name)}\n\n")
            for item in items:
                if i < 3:
                    f.write(f".. autoclass:: {mod}.{item}\n\n")
                else:
                    f.write(f".. autofunction:: {mod}.{item}\n\n")


# Auto-generate isothermal equation of state documentation with available strings
output_rst = os.path.join(
    os.path.dirname(__file__), "autogenerated/isothermal_eos.rst"
)

with open(output_rst, "w") as f:

    for mod, clses in isothermal_equations_of_state_dict.items():

        module = __import__(mod, fromlist=["dummy"])
        for cls in clses:

            # Skip base classes
            if cls[-4:].lower() == "base":
                continue

            obj = getattr(module, cls)

            string_key = None
            if obj in eos_methods.values():
                for key, value in eos_methods.items():
                    if value == obj:
                        string_key = key
                        break

            if string_key is not None:
                cls_name = eos_names[string_key]
                heading = f"{cls_name} (\"{string_key}\")"
                f.write(f"{heading}\n{'^' * len(heading)}\n\n")
            else:
                f.write(f"{cls}\n{'^' * len(cls)}\n\n")

            f.write(f".. autoclass:: {mod}.{cls}\n")
            f.write("   :no-members:\n")
            f.write("   :no-undoc-members:\n")
            f.write("   :no-inherited-members:\n")
            f.write("   :no-index:\n\n")

# Auto-generate thermal equation of state documentation with available strings
output_rst = os.path.join(
    os.path.dirname(__file__), "autogenerated/thermal_eos.rst"
)

with open(output_rst, "w") as f:

    for mod, clses in thermal_equations_of_state_dict.items():

        module = __import__(mod, fromlist=["dummy"])
        for cls in clses:

            # Skip base classes
            if cls[-4:].lower() == "base":
                continue

            obj = getattr(module, cls)

            string_key = None
            if obj in eos_methods.values():
                for key, value in eos_methods.items():
                    if value == obj:
                        string_key = key
                        break

            if string_key is not None:
                cls_name = eos_names[string_key]
                heading = f"{cls_name} (\"{string_key}\")"
                f.write(f"{heading}\n{'^' * len(heading)}\n\n")
            else:
                f.write(f"{cls}\n{'^' * len(cls)}\n\n")

            f.write(f".. autoclass:: {mod}.{cls}\n")
            f.write("   :no-members:\n")
            f.write("   :no-undoc-members:\n")
            f.write("   :no-inherited-members:\n")
            f.write("   :no-index:\n\n")

# Auto-generate property modifier documentation
output_rst = os.path.join(
    os.path.dirname(__file__), "autogenerated/property_modifiers.rst"
)

with open(output_rst, "w") as f:
    for string, func in modifier_functions.items():
        name = modifier_names[string]
        header = f"{name} (available as \"{string}\")"
        f.write(f"{header}\n")
        f.write('"' * len(header) + "\n\n")
        f.write(f".. autofunction:: burnman.eos.property_modifiers.{func.__name__}\n")
        f.write("   :no-index:\n\n")


# Auto-generate solution model documentation
output_rst = os.path.join(
    os.path.dirname(__file__), "autogenerated/solution_models.rst"
)

# All classes are contained in burnman.classes.solutionmodel
module_name = burnman.classes.solutionmodel.__name__
classes = [
    (name, obj) for name, obj in burnman.classes.solutionmodel.__dict__.items()
    if inspect.isclass(obj) and obj.__module__ == module_name
]

with open(output_rst, "w") as f:
    for name, obj in classes:
        # Only consider classes
        if isinstance(obj, type):

            # Skip base classes
            if name.endswith("Base") or name.endswith("Model") or name == "Interaction":
                continue

            f.write(f"{name}\n" + '"' * len(header) + "\n\n")
            f.write(f".. autoclass:: {module_name}.{name}\n")
            f.write("   :no-members:\n")
            f.write("   :no-undoc-members:\n")
            f.write("   :no-inherited-members:\n")
            f.write("   :no-index:\n\n")

# Auto-generate elastic solution model documentation
output_rst = os.path.join(
    os.path.dirname(__file__), "autogenerated/elastic_solution_models.rst"
)

module_name = burnman.classes.elasticsolutionmodel.__name__
classes = [
    (name, obj) for name, obj in burnman.classes.elasticsolutionmodel.__dict__.items()
    if inspect.isclass(obj) and obj.__module__ == module_name
]

with open(output_rst, "w") as f:
    for name, obj in classes:
        # Only consider classes
        if isinstance(obj, type):

            # Skip base classes
            if name.endswith("Base") or name.endswith("Model") or name == "Interaction":
                continue

            f.write(f"{name}\n" + '"' * len(name) + "\n\n")
            f.write(f".. autoclass:: {module_name}.{name}\n")
            f.write("   :no-members:\n")
            f.write("   :no-undoc-members:\n")
            f.write("   :no-inherited-members:\n")
            f.write("   :no-index:\n\n")

# Auto-generate averaging scheme documentation
output_rst = os.path.join(
    os.path.dirname(__file__), "autogenerated/averaging_schemes.rst"
)

module_name = burnman.classes.averaging_schemes.__name__
classes = [
    (name, obj) for name, obj in burnman.classes.averaging_schemes.__dict__.items()
    if inspect.isclass(obj) and obj.__module__ == module_name
]

with open(output_rst, "w") as f:
    for name, obj in classes:
        # Only consider classes
        if isinstance(obj, type):

            # Skip base class
            if name == "AveragingScheme":
                continue

            f.write(f"{name}\n" + '"' * len(name) + "\n\n")
            f.write(f".. autoclass:: {module_name}.{name}\n")
            f.write("   :no-members:\n")
            f.write("   :no-undoc-members:\n")
            f.write("   :no-inherited-members:\n")
            f.write("   :no-index:\n\n")

Overview
========

Overall Structure
-----------------

BurnMan is designed to be a general mineral physics and seismological toolkit
which can enable a user to calculate (or fit) the physical and chemical properties
of endmember minerals, fluids/melts, solid solutions, and composite assemblages.
Such properties include:

* the thermodynamic free energies, allowing phase equilibrium calculations,
  endmember activities, chemical potentials and oxygen (and other) fugacities.
* entropy, enabling the user to calculate isentropes for a given assemblage.
* volume, to allow the user to create density profiles.
* seismic velocities, including Voigt-Reuss-Hill and Hashin-Strikman bounds
  and averages.

Data and functions are provided to allow the user to compare calculated
isentropes and seismic velocity profiles  to profiles computed for other
compositions or constrained by seismology.

BurnMan is written in the Python language and is run from the command
line.  This allows the library to be incorporated into other projects.
BurnMan makes extensive use of `SciPy <http://www.scipy.org/>`_ and `NumPy <http://www.numpy.org/>`_, which are widely used Python
libraries for scientific computation.  `Matplotlib <http://matplotlib.org/>`_ is used to display results
and produce publication quality figures.  The computations are consistently
formulated in terms of SI units.

The toolkit includes:
 
* the full codebase, including many equations of state and solution models
* popular datasets already coded into burnman-usable format
* a tutorial on the basic use of BurnMan
* a large collection of annotated examples
* an extensive suite of unit tests to ensure code functions as intended
* a series of benchmarks comparing BurnMan output with published data
* a directory containing user-contributed code from published papers
 
This software has been designed to allow the end-user a great deal of freedom
to do whatever calculations they may wish. We have endeavoured to provide
examples and benchmarks which cover the most popular uses of the software,
some of which are included in the figure below. This list is certainly not
exhaustive, and we will definitely have missed interesting
applications. As a result we will be very happy to accept contributions in
form of corrections, examples, or new features.

.. graphviz::

   digraph workflows {
     rankdir=TB;
     experimental_data [shape=ellipse,label="experimental data"]
     assemblage [shape=ellipse,label="assemblage"]
     bulk [shape=ellipse,label="bulk composition"]
     pt [shape=ellipse,label="P-T conditions"]
     mineral [shape=diamond,label="mineral",fillcolor=red, style=filled]
     geotherm [shape=diamond,label="geotherm",fillcolor=red,
     style=filled]
     seismic_model [shape=diamond,label="seismic.SeismicTable",fillcolor=red, style=filled]
     solution [shape=diamond,label="solution",fillcolor=red, style=filled]
     composite[shape=diamond,label="composite",fillcolor=red,
     style=filled]
     av
     [shape=diamond,label="averaging_schemes.AveragingScheme",fillcolor=red,
     style=filled]
     fit [shape=box,label="tools.fit_PVT_data",fillcolor=yellow, style=filled]
     eqm
     [shape=box,label="equilibriumassemblage.gibbs_minimizer",fillcolor=yellow,
     style=filled]
     chempot [shape=box,label="chemicalpotentials.chemical_potentials",fillcolor=yellow,
     style=filled]
     subgraph experimental {
         label="Equation of state fitting";
         mineral ->fit;
         experimental_data ->fit;
         fit-> "optimized EoS parameters";
      }
      subgraph thermodynamic {
         label="Thermodynamic calculations";
         mineral -> composite;
         solution -> composite;
         bulk -> eqm;
	 pt -> eqm;
         composite -> eqm;
         eqm -> assemblage;
         assemblage -> chempot;
	 pt -> chempot;
	 chempot -> "chemical potentials / fugacities";
      }
      subgraph seismic {
         label="Seismic comparisons";
	 geotherm -> pt;
	 assemblage ->"seismic velocities";
         av -> "seismic velocities";
	 "seismic velocities" -> "velocity comparison"
	 seismic_model-> "velocity comparison"
      }
   }

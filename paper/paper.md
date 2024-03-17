---
title: BurnMan -- a Python toolkit for planetary geophysics, geochemistry and thermodynamics
tags:
  - Python
  - planetary
  - geophysics
  - geochemistry
  - thermodynamics
authors:
  - name: Robert Myhill
    orcid: 0000-0001-9489-5236
    corresponding: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Sanne Cottaar
    orcid: 0000-0003-0493-6570
    affiliation: 2
  - name: Timo Heister
    orcid: 0000-0002-8137-3903
    affiliation: 3
  - name: Ian Rose
    affiliation: 4
  - name: Cayman Unterborn
    orcid: 0000-0001-8991-3110
    affiliation: 5
  - name: Juliane Dannberg
    orcid: 0000-0003-0357-7115
    affiliation: 6
  - name: Rene Gassmoeller
    orcid: 0000-0001-7098-8198
    affiliation: 6
affiliations:
 - name: University of Bristol, UK
   index: 1
 - name: University of Cambridge, UK
   index: 2
 - name: Clemson University, USA
   index: 3
 - name: Independent Researcher, USA
   index: 4
 - name: Southwest Research Institute, USA
   index: 5
 - name: University of Florida, USA
   index: 6
date: 20 September 2022
bibliography: paper.bib


# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
Many branches of physics, chemistry and Earth sciences build
complex materials from simpler constructs: rocks are composites
made up of one or more phases; phases are solutions of more than one
endmember; and endmembers are usually mixtures of more than one element.
The properties of the endmember building blocks at different pressures and
temperatures can be modelled using a wide array of different equations of state.
There are also many models for the averaging of endmember properties
within solutions and composite materials. Once calculated,
the physical properties of composite materials can be used in many
different ways.

`BurnMan` is an open source, extensible mineral physics module
written in Python. It implements several different methods to
calculate the physical properties of natural materials.
The toolbox has a class-based, modular
design that allows users to calculate many low-level properties that
are not easily accessed using existing codes, and to combine various tools in
novel, creative ways. The module includes:

* over a dozen static and thermal equations of state for pure phases;
* commonly-used solution model formalisms (ideal, (a)symmetric, subregular)
and a formalism that allows users to define their own excess energy functions;
* popular endmember and solution datasets for solids and melts,
including @Holland:2011, @deKoker:2013 and @Stixrude:2022;
* an anisotropic equation of state [@Myhill:2022];
* a consistent method for combining phases into a composite assemblage,
with seismic averaging schemes including Voigt, Reuss, Voigt-Reuss-Hill
and the Hashin-Shtrikman bounds;
* a common set of methods to output thermodynamic and thermoelastic
properties for all materials;
* a solver to chemically equilibrate composite materials;
* optimal least squares fitting routines for multivariate experimental
data with (potentially correlated) errors. These allow (for example) simultaneous
fitting of pure phase and solution model parameters to experimental volumes,
seismic velocities and enthalpies of formation;
* "Planet" and "Layer" classes that self-consistently calculate
gravity, pressure, density, mass, moment of inertia and seismic velocity profiles
given chemical, thermal and dynamic constraints;
* geothermal profiles from the literature as well as the option to calculate
adiabatic profiles based on mineral assemblage;
* a set of high-level functions which create files readable by
seismological and geodynamic software, including: Mineos [@Masters:2011],
AxiSEM [@NissenMeyer:2014] and
ASPECT [@Kronbichler:2012;@aspect-doi-v2.4.0;@aspectmanual]; and
* a `Composition` class, which provides a framework to convert between mass, molar,
and elemental compositions, convert to different chemical component systems,
and add or subtract components. 

The project includes over 40 annotated examples,
an extensive suite of unit tests and benchmarks, and 
a directory containing user-contributed code from published papers.
A multipart tutorial illustrates key functionality, including the
functions required to create the figures in this paper
(\url{https://burnman.readthedocs.io/en/latest/tutorial.html}).
Using `BurnMan` requires only moderate Python skills, and its modular nature
means that it can easily be customised. 

# Statement of need
Earth, planetary and materials scientists are interested in a number
of different material properties, including seismic velocities,
densities and heat capacities as functions of pressure and temperature.
Many of these properties are connected to each
other by physical laws such as Maxwell's relations.
Building models of individual phases to compute
these properties can be time-consuming and prone to error.
It is desirable to have well-tested and benchmarked software that
makes it easy to calculate the properties of complex
composite materials from existing models, and to parameterize
new models from experimental data.
Furthermore, there are many common scientific workflows that require
thermodynamic and thermoelastic properties as input.
These are the needs satisfied by the `BurnMan` module.

# The BurnMan project
The focus of `BurnMan` was originally on the seismic properties of the
lower mantle [@Cottaar:2014]. Its scope has now expanded
to encompass the thermodynamic and thermoelastic properties of any
geological and planetary materials
(see \url{https://github.com/geodynamics/burnman/releases}
for the history of improvements). Pure phase equations of state are designed
to be sufficiently flexible to model real-world materials from the Earth's core
to the shallowest parts of the crust (e.g. \autoref{fig:qtzproperties}).
Solution model formulations with varying levels
of non-ideality are included (e.g. \autoref{fig:garnetsolution}),
including both Gibbs and Helmholtz formulations [@Myhill:2018].
Functions are provided to convert solution models from one
endmember basis to another [@Myhill:2021]. A `Composite` class
allows calculation of the properties of assemblages
containing several phases and includes several seismic averaging schemes.

![Heat capacity and bulk sound velocities of quartz through the alpha-beta
quartz transition as found in [@Stixrude:2011]. This transition is modelled
via a Landau-type model. 
\label{fig:qtzproperties}](figures/quartz_properties.png)

![Properties of pyrope-grossular garnet at 1 GPa according to a published
model [@Jennings:2015], as output by `BurnMan`.
The excess Gibbs energy is useful for calculating phase equilibria by
Gibbs minimization, while the endmember activities can be used to
determine equilibrium via the equilibrium relations [@Holland:1998].
\label{fig:garnetsolution}](figures/mg_ca_gt_properties.png)

`BurnMan` also includes planetary `Layer` and `Planet` classes
that can be used to construct planetary models with self-consistent
pressure, gravity and density profiles
and calculate seismic properties through those bodies.
\autoref{fig:zog} shows the output from a model that resembles Earth.
Tools are provided to compare predicted seismic properties
with published seismic models of the Earth,
and to produce input files to compute synthetic
seismic data using other codes, including
AxiSEM [@NissenMeyer:2014] and Mineos [@Woodhouse:1988;@Masters:2011].

![A 1D profile through Planet Zog, a planet much like Earth, with an
inner and outer core (blue and orange layers),
isentropic convecting lower and upper mantle (green and red), and a
depleted lithosphere (lilac) split into mantle and crust.
The mineralogy/composition of each layer is chosen by the user.
Zog has the same mass (5.972e+24 kg) and moment of inertia factor
(0.3307) as Earth. `BurnMan` ensures that the gravity and pressure
profiles satisfy hydrostatic equilibrium, and allows different layers
to have different thermal profiles, including an isentropic profile
with thermal boundary layers
(shown here for the upper mantle, lower mantle and for the core).
Depth dependent changes to density, gravity, pressure (solid blue lines)
are compared with the Preliminary Reference Earth Model
(PREM; dotted orange line, [@Dziewonski:1981]). 
The computed geotherm is compared to several from the literature
[@Stacey:1977;@Brown:1981;@Anderson:1982;@Alfe:2007;@Anzellini:2013].
\label{fig:zog}](figures/zog.png)


`BurnMan` also includes many utility functions. These include functions
that fit parameters for pure phase
and solution models to experimental data
including full error propagation (\autoref{fig:fit}). 
Other fitting functions include `fit_composition_to_solution()` and  
`fit_phase_proportions_to_bulk_composition()` that
use weighted constrained least squares using cvxpy [@Diamond:2016] to
estimate endmember or phase proportions given a bulk composition.
These fitting functions apply appropriate non-negativity constraints 
(i.e. that no species can have negative proportions on a site, 
and that no phase can have a negative abundance
in the bulk). An example of `fit_phase_proportions_to_bulk_composition()`
that uses real experimental data [@Bertka:1997] is shown in \autoref{fig:mars_fit}.
Loss of an independent component from the bulk composition
can be tested by adding another phase with the composition of
that component (e.g. Fe) and checking that the amount of that phase
is zero within uncertainties. 

![Optimized fit of a PV equation of state [@Holland:2011] to
stishovite data [@Andrault:2003], including 95% confidence intervals
on both the volume and the bulk modulus.
\label{fig:fit}](figures/stishovite_fitting.png)

![Mineral phase proportions in the mantle of Mars, estimated by using 
the method of constrained least squares on high pressure experimental
data [@Bertka:1997]. Weighted residuals
(misfits) are also shown, indicating that the majority of experimental
run products are consistent with the reported starting composition.
\label{fig:mars_fit}](figures/fit_mars_experiments.png)

`BurnMan` does not attempt to replicate Gibbs minimization codes,
of which there are many, such as PerpleX [@Connolly:2009],
MELTS [@Ghiorso:1995],
MageMin [@Riel:2022],
TheriakDomino [@deCapitani:2010],
HeFESTo [@Stixrude:2022] and
FactSAGE [@Bale:2002].
Instead, it provides two methods to deal with the problem of
thermodynamic equilibrium: (1) reading in a pressure-temperature table
of precalculated properties into a Material class that allows derivative
properties to be calculated, and (2) a function called `equilibrate()`
that equilibrates a known assemblage under user-defined constraints.
This function requires an assemblage
(e.g. olivine, garnet and orthopyroxene), a starting bulk composition,
desired equality constraints, and optionally one or more
compositional degrees of freedom. The equilibrate function solves the
equilibrium relations [@Holland:1998] using a
damped Newton root finder [@Nowak:1991].

The `equilibrate()` function allows the user to select from a number of equality
constraints, including fixed pressure, temperature, entropy or volume, or
compositional equalities such as a fixed molar fraction of a phase,
or a certain ratio of Mg and Fe on a particular site. 
The number of constraints required is two at fixed bulk composition, 
and one more for each degree of compositional freedom. An example of the use
of the equilibrate function is shown in \autoref{fig:eqm}.
Full details may be found in the manual and tutorial.

![The olivine phase diagram at three different temperatures as computed
using the equilibrate routines in `BurnMan`. The solution model properties
are taken from the literature [@Stixrude:2011].
\label{fig:eqm}](figures/olivine_equilibrium.png)

# Past and ongoing research projects
In addition to mantle studies
[@Cottaar:2014b;@Ballmer:2017;@Jenkins:2017;@Thomson:2019;@Houser:2020],
`BurnMan` has been used to investigate Earth's core [@Irving:2018],
in phase equilibria studies [@Myhill:2017;@Ishii:2019], to develop new models
for anisotropic thermodynamics [@Myhill:2022], to constrain the interiors
of exoplanets [@Unterborn:2016;@Unterborn:2019],
and to provide input for geodynamic simulations
[@Heister:2017;@Dannberg:2021].

# Acknowledgments
BurnMan was initiated at, and follow-up research support was received through,
CIDER (NSF FESD grant 1135452). The authors have been supported by the
Computational Infrastructure for Geodynamics initiative (CIG),
through the National Science Foundation (U.S.) under Award No. EAR-0949446.
They have also received support from The University of California-Davis.

Robert Myhill was supported by the Science and Technologies Funding Council
(U.K.) under Award No. ST/R001332/1 and through the
Natural Environment Research Council (U.K.) Large Grant MC-squared
(Award No. NE/T012633/1). He would also like to thank M. Ghiorso for useful
discussions and B. Watterson for the concept of Planet Zog.

Sanne Cottaar has received funding from the European Research Council (ERC)
under the European Union's Horizon 2020 research and innovation programme
(ZoomDeep; Award No. 804071). This funding covered collaborative visits.

Timo Heister was partially supported by NSF Award DMS-2028346,
OAC-2015848, EAR-1925575, and by the Computational Infrastructure in
Geodynamics initiative (CIG), through the NSF under Award EAR-0949446
and EAR-1550901 and The University of California -- Davis.

Rene Gassmoeller and Juliane Dannberg were supported by NSF Awards
EAR-1925677 and EAR-2054605, and by the Computational Infrastructure for
Geodynamics (CIG) through the NSF under Award EAR-0949446 and EAR-1550901
and the University of California – Davis.

The BurnMan code has been contributed to the Computational Infrastructure
for Geodynamics (CIG) and is hosted at geodynamics.org.

# References

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
 - name: University of Utah, USA
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
endmember and endmembers are usually mixtures of more than one element.
The properties of the endmember building blocks at different pressures and
temperatures can be provided by a wide array of different equations of state.
The averaging of the endmember properties within solutions and composite
materials can also be achieved in several different ways.

`BurnMan` is an open source, extensible mineral physics toolbox designed to
calculate the physical properties of natural materials
within the Earth and other planets. The toolbox has a class-based, modular
design that allows users to calculate many low-level properties that
are not accessible using existing codes, and to combine various tools in
novel, creative ways. The module includes:

* many static and thermal equations of state for endmembers
* several solution model formalisms (ideal, (a)symmetric, subregular)
and an option for user-defined functional forms
* popular endmember and solution datasets 
(including [@Holland:2011] and [@Stixrude:2011])
* a consistent method for combining phases into a composite assemblage,
with several seismic averaging schemes
* a common set of methods to output thermodynamic and thermoelastic
properties for all materials
* an equilibrium reaction solver for composites
* optimal least squares fitting routines for multivariate experimental
data with (potentially correlated) errors. As an example, such functions
can be used to simultaneously fit volumes, seismic velocities
and enthalpies
* a "Planet" class, which self-consistently calculates gravity profiles,
mass, moment of inertia and temperature structure of planets given
appropriate chemical and dynamic constraints
* published geotherms
* a set of high-level functions which create files readable by
seismological and geodynamic software, including: Mineos [@Masters:2011],
AxiSEM [@NissenMeyer:2014] and
ASPECT [@Kronbichler:2012;@aspect-doi-v2.4.0;@aspectmanual]

The project also includes a tutorial, a large collection of annotated examples,
an extensive suite of unit tests and benchmarks, and 
a directory containing user-contributed code from published papers.
Use of the code requires only moderate Python skills, and its modular nature
means that it can easily be customised.
The `BurnMan` project can be found at
\url{https://github.com/geodynamics/burnman}.

# Statement of need
Earth Scientists are interested in a number of different material properties,
including seismic velocities, heat capacities and densities as functions of
pressure and temperature. Many of these properties are connected to each
other by physical laws. Building models of individual phases to compute
these properties and checking them can be time-consuming and prone to error,
so it is desirable to have well-tested and benchmarked software that
provides convenient functions to calculate the properties of complex
composite materials, and to fit new models.
Furthermore, there are many common workflows that require physical properties
as input, or refine material parameters to fit experimentally-determined
physical properties. These are the needs satisfied by the `BurnMan` module.

# The BurnMan project
When `BurnMan` was first released [@Cottaar:2014], its focus was on
the seismic properties of the lower mantle, using a single endmember
mineral database [@Stixrude:2011] as a foundation. Since then,
its scope has expanded considerably. `BurnMan` now contains equations
of state for minerals and melts from several published datasets.
A common set of methods for all equations of state allows easy
access to many thermodynamic properties
(e.g. \autoref{fig:qtzproperties}).

![Heat capacity and bulk sound velocities of quartz through the alpha-beta
quartz transition as found in [@Stixrude:2011]. This transition is modelled
via a Landau-type model.
\label{fig:qtzproperties}](figures/quartz_properties.png)

An extension to anisotropic properties has recently been added
[@Myhill:2022]. Several solution model formulations are included
(with one example model shown in \autoref{fig:garnetsolution}),
including methods to convert from one endmember basis to another
[@Myhill:2021]. A composite material model
allows calculation of the properties of assemblages of several phases
and includes several seismic averaging schemes.

![Properties of pyrope-grossular garnet at 1 GPa according to a published
model [@Jennings:2015], as output by `BurnMan`.
The excess Gibbs energy is useful for calculating phase equilibria by
Gibbs minimization, while the endmember activities can be used to
determine equilibrium via the equilibrium relations [@Holland:1998].
\label{fig:garnetsolution}](figures/mg_ca_gt_properties.png)

`BurnMan` includes many higher level functions that build on
its material model classes. These include functions to
fit thermodynamic models for both endmembers
and solutions including full error propagation (\autoref{fig:fit}),
construction of self-consistent planetary models (\autoref{fig:zog})
and the calculation of seismic properties. Tools are provided to
compare the seismic predictions with published seismic models of
the Earth, and to produce input files to compute synthetic
seismic data using other codes, including
AxiSEM [@NissenMeyer:2014] and Mineos [@Woodhouse:1988;@Masters:2011].

![Optimized fit of a PV equation of state [@Holland:2011] to
stishovite data [@Andrault:2003], including 95% confidence intervals
on both the volume and the bulk modulus.
\label{fig:fit}](figures/stishovite_fitting.png)

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
The computed geotherm is compared to several from the literature
[@Stacey:1977;@Brown:1981;@Anderson:1982;@Alfe:2007;@Anzellini:2013].
The other properties are compared to the PREM [@Dziewonski:1981].
\label{fig:zog}](figures/zog.png)

`BurnMan` also provides many utility classes and functions, including
those required to undertake bulk and phase chemical calculations,
to aid preparation and analysis of laboratory experiments.
An example that uses bulk and phase compositions from samples
to estimate phase proportions is shown in \autoref{fig:mars_fit}.

![Mineral phase proportions in the mantle of Mars, estimated from
high pressure experimental data [@Bertka:1997]. Weighted residuals
(misfits) are also shown, indicating that the majority of experimental
run products are consistent with the reported starting composition.
\label{fig:mars_fit}](figures/fit_mars_experiments.png)

`BurnMan` does not attempt to replicate Gibbs minimization codes,
of which there are many, such as PerpleX [@Connolly:2009],
MageMin [@Riel:2022],
TheriakDomino [@deCapitani:2010],
HeFESTo [@Stixrude:2022] and
FactSAGE [@Bale:2002].
Instead, it provides two methods to deal with the problem of
thermodynamic equilibrium: (1) reading in a pressure-temperature table
of precalculated properties into a Material class
(from which derivative properties can be calculated),
and (2) an equilibrate function
that chemically equilibrates a known assemblage under constraints
(two or three choices from fixed pressure, temperature, entropy, volume,
phase proportions and compositions).
An example is shown in \autoref{fig:eqm}.
The equilibrate function is similar to that implemented in ThermoCalc
[@Holland:1998].

![The olivine phase diagram at three different temperatures as computed
using the equilibrate routines in `BurnMan`. The solution model properties
are taken from the literature [@Stixrude:2011].
\label{fig:eqm}](figures/olivine_equilibrium.png)

The features used in `BurnMan` are documented in the codebase
(\url{https://github.com/geodynamics/burnman}), and in the
manual (\url{https://burnman.readthedocs.io}). A Jupyter notebook tutorial
(available at \url{https://geodynamics.github.io/burnman/#download})
introduces some of the core functionality, including demonstrations that
create the figures presented here. More features are investigated
in the extensive set of examples.

# Past and ongoing research projects
In addition to mantle studies
[@Cottaar:2014b;@Ballmer:2017;@Jenkins:2017;@Thomson:2019;@Houser:2020],
`BurnMan` has been used to investigate Earth's core [@Irving:2018],
in phase equilibria studies [@Myhill:2017;@Ishii:2019], to develop new models
for anisotropic thermodynamics [@Myhill:2022], to constrain the interiors
of exoplanets [@Unterborn:2016;@Unterborn:2019],
and to provide input for geodynamic simulations
[@Heister:2017;@Dannberg:2021].

# Acknowledgements
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

Timo Heister was partially supported by NSF Award DMS-1522191, DMS-1820958,
OAC-1835452, by the Computational Infrastructure in Geodynamics initiative
(CIG), through the NSF under Award EAR- 0949446 and EAR-1550901 and
The University of California – Davis, and by Technical Data Analysis, Inc.
through US Navy SBIR N16A-T003.

Rene Gassmoeller and Juliane Dannberg were supported by NSF Awards
EAR-1925677 and EAR-2054605, and by the Computational Infrastructure for
Geodynamics (CIG) through the NSF under Award EAR-0949446 and EAR-1550901
and the University of California – Davis.

The BurnMan code has been contributed to the Computational Infrastructure
for Geodynamics (CIG) and is hosted at geodynamics.org.

# References
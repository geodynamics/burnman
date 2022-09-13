---
title: BurnMan -- a Python toolkit for planetary geophysics, geochemistry and thermodynamics'
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
date: 13 September 2022
bibliography: paper.bib


# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
Many branches of physics, chemistry and Earth sciences conceptually build
complex materials up from simpler constructs: rocks are composites
made up of one or more phases; phases are solutions of more than one endmember
and endmembers are usually mixtures of more than one atom / oxide. The
properties of the endmember building blocks at different pressures and
temperatures can be provided by a wide array of different equations of state.
The averaging of those properties within solutions and composite materials
can also be achieved in several different ways.

# Statement of need
Earth Scientists are interested in a number of different material properties,
including seismic velocities, heat capacities and density as a function of
pressure and temperature. Many of these properties are connected to each
other by physical laws. Building up models of
individual phases to compute these properties and checking them can be time-consuming and prone to error,
and so it is desirable to have well-tested and benchmarked software that
provides convenient functions to calculate the properties of complex
composite materials, and to fit new models.

`BurnMan` is an open source, extensible mineral physics toolbox designed to
calculate the physical properties of natural materials
within the Earth and other planets. The toolbox has a class-based, modular
design that allows users to calculate many low-level properties that
are not accessible using existing codes.

When `BurnMan` was first released [@Cottaar:2014], its focus was on
the seismic properties of the lower mantle, using a single endmember mineral
database [@Stixrude:2011] as a foundation. Since then,
its scope has expanded considerably. Some workflows of the current code
are provided in \autoref{fig:workflows}. 

![Some possible BurnMan workflows. Input/output shown as ellipses, functions as
rectangles and classes as diamonds.\label{fig:workflows}](figures/workflows.png)

`BurnMan` now contains equations of state
for minerals and melts from several published datasets. A model for anisotropy
has recently been added. Several solution model formulations are included
(\autoref{fig:garnetsolution}),
and a composite material model that allows several seismic averaging schemes.


![Properties of pyrope-grossular garnet at 1 GPa according to a published
model [@Jennings:2015], as output by `BurnMan`.
\label{fig:garnetsolution}](figures/mg_ca_gt_properties.png)


These material models can be used for many higher levels purposes,
including the fitting of thermodynamic models for both endmembers
and solutions including full error propagation (\autoref{fig:fit}),
construction of self-consistent planetary models
and the calculation of seismic properties. An example of the
planet builder in `BurnMan` is shown in \autoref{fig:zog}.

![Optimized fit of a PVT equation of state to periclase data, including
95% confidence intervals.\label{fig:fit}](figures/example_fit_eos15.png)

![A 1D profile through Planet Zog, a planet much like Earth, with an
inner and outer core, convecting mantle (with thermal boundary layers),
depleted lithospheric mantle and crust. Zog has the same mass (5.972e+24)
and moment of inertia factor (0.3307) as Earth. Properties calculated
self-consistently using BurnMan. The computed geotherm is compared to
several from the literature
[@Stacey:1977;@Brown:1981;@Anderson:1982;@Alfe:2007;@Anzellini:2013].
\label{fig:zog}](figures/zog.png)

Utility classes and functions are also provided to undertake bulk and 
phase chemical calculations, to aid preparation and analysis of experiments.

`BurnMan` does not attempt to replicate Gibbs minimization codes,
of which there are many (PerpleX, HeFESTo, TheriakDomino, FactSAGE and MageMin).
Instead, it provides two methods to deal with the problem of
thermodynamic equilibrium: (1) reading in a P-T table of precalculated
properties into a Material class, and (2) a function (`burnman.equilibrate`)
that chemically equilibrates a known assemblage under constraints
(two or three choices from fixed pressure, temperature, entropy, volume,
phase proportions and compositions).
An example is shown in Figure \autoref{fig:eqm}.
The equilibrate function is similar to that implemented in ThermoCalc
[@Holland:1998].

![The olivine phase diagram at three different temperatures.
\label{fig:eqm}](figures/olivine_equilibrium.png)

The features used in `BurnMan` are documented in the codebase
(\url{https://github.com/geodynamics/burnman}), and in the
manual (\url{https://burnman.readthedocs.io}). A Jupyter notebook tutorial
(available at \url{https://geodynamics.github.io/burnman/#download})
introduces some of the core functionality, and more features are investigated
in the extensive set of examples.

# Past and ongoing research projects
In addition to lower mantle studies [@Cottaar:2014b;@Jenkins:2017],
`BurnMan` has been used to investigate Earth's core [@Irving:2018],
in phase equilibria studies [@Myhill:2017;@Ishii:2019], to develop new models
for anisotropic thermodynamics [@Myhill:2022], to constrain the interiors
of exoplanets [@Unterborn:2016], and to provide input for geodynamic simulations
[@Heister:2017;@Dannberg:2021]. 

# Acknowledgements
We thank all the members of the CIDER Mg/Si team for their input:
Valentina Magni, Yu Huang, JiaChao Liu, Marc Hirschmann,
and Barbara Romanowicz. We also thank Lars Stixrude for providing benchmarking
calculations and Zack Geballe, Motohiko Murakami, Bill McDonough,
Quentin Williams, Wendy Panero, and Wolfgang Bangerth for helpful discussions.

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

Timo Heister was partially supported by NSF Award DMS-1522191, DMS-1820958,
OAC-1835452, by the Computational Infrastructure in Geodynamics initiative
(CIG), through the NSF under Award EAR- 0949446 and EAR-1550901 and
The University of California – Davis, and by Technical Data Analysis, Inc.
through US Navy SBIR N16A-T003.

The BurnMan code has been contributed to the Computational Infrastructure
for Geodynamics (CIG) and is hosted at geodynamics.org.

# References

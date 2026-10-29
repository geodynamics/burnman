.. _ref-seismic-models:

Seismic Models
==============

Overview
--------

BurnMan allows for direct visual and quantitative comparison with seismic velocity models.
Various ways of plotting can be found in the examples.
Quantitative misfits between two profiles include an L2-norm and a chi-squared misfit, but user defined norms can be implemented.
A seismic model in BurnMan is
an object that provides pressure, density, and seismic velocities (:math:`V_P, V_\Phi, V_S`) as a function of depth.

See :ref:`ref-api-seismic` for a list of built-in seismic models available in BurnMan.

Comparison to Seismic profiles
------------------------------

To compare to seismically constrained profiles, BurnMan provides the 1D seismic velocity model PREM :cite:`dziewonski1981`.
One can choose to evaluate :math:`V_P, V_\Phi, V_S, \rho, K_S` and/or :math:`G`.
The user can input their own seismic profile, an example of which is included using AK135 :cite:`kennett1995`.

Besides standardized 1D radial profiles, one can also compare to regionalized average profiles for the lower mantle.
This option accommodates the observation that the lowermost mantle can be clustered into two regions, a 'slow' region, which represents the so-called Large Low Shear Velocity Provinces, and 'fast' region, the continuous surrounding region where slabs might subduct :cite:`Lekic2012`.
This clustering as well as the averaging of the 1D model occurs over five tomographic S wave velocity  models (SAW24B16: :cite:`megnin2000`; HMSL-S: :cite:`houser2008`; S362ANI: :cite:`kustowski2008`; GyPSuM: :cite:`Simmons2010`; S40RTS: :cite:`Ritsema2011`).
The strongest deviations from PREM occur in the lowermost 1000 km.
Using the 'fast' and 'slow' S wave velocity profiles is therefore most important when interpreting the lowermost mantle. Suggestion of compositional variation between these regions comes from seismology :cite:`to2005,He2012` as well as geochemistry :cite:`Deschamps2012,jackson2010`.
Based on thermo-chemical convection models, :cite:`Styles2011` also show that averaging profiles in thermal boundary layers may cause problems for seismic interpretation.

We additionally apply cluster analysis to and provide models for P wave velocity based on two tomographic models (MIT-P08: :cite:`Li2008`; GyPSuM: :cite:`Simmons2012`).
The clustering results correlate well with the fast and slow regions for S wave velocities; this could well be due to the fact that the initial model for the P wave velocity models is scaled from S wave tomographic velocity models.
Additionally, the variations in P wave velocities are a lot smaller than for S waves.
For this reason using these adapted models is most important when comparing the S wave velocities.

While interpreting lateral variations of seismic velocity in terms of composition and temperature is a major goal :cite:`Trampert2004,Mosca2012`, to determine the bulk composition the current challenge appears to be concurrently fitting absolute P and S wave velocities and incorporate the significant uncertainties in mineral physical parameters).

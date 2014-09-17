

Methods
=======

.. _sec:EoS:
Calculating Thermoelastic Properties
------------------------------------


To calculate the bulk (:math:`K`) modulus, shear modulus (:math:`G`) and density (:math:`\rho`) of a material at a given pressure (:math:`P`) and temperature (:math:`T`), optionally defined by a geotherm) and determine the seismic velocities (:math:`V_S, V_P, V_\Phi`), one uses an Equation of State (EoS).
Currently the following EoSs are supported in BurnMan: the Birch-Murnaghan formulation (excludes temperature effects) :cite:`Poirier1991`, and the Birch-Murnaghan formulation with a Mie-Grüneisen-Debye temperature correction as formulated by :cite:`Stixrude2005`.
To calculate these thermoelastic parameters, the EoS requires the user to input three parameters: pressure, temperature, the phases and their molar fractions.
These inputs and outputs are further discussed in Section :ref:`input`.



Isothermal calculations: Birch-Murnaghan
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Birch-Murnaghan equation is an isothermal Eulerian finite-strain EoS relating pressure and volume.
The negative finite-strain (or compression) is defined as


.. math::
    f=\frac{1}{2}\left[\left(\frac{V}{V_0}\right)^{-2/3}-1\right],
    :label: f

where :math:`V` is the volume at a given pressure and :math:`V_0` is the volume at a reference state (:math:`P = 10^5` Pa , :math:`T` = 300 K).
The pressure and elastic moduli are derived from a third order Taylor expansion of Helmholtz free energy in :math:`f` and evaluating the appropriate volume and strain derivatives (see, e.g., :cite:`poirier1991`).
For an isotropic material one obtains for the pressure, isothermal bulk modulus, and shear modulus:


.. math::
   P=3K_0f\left(1+2f\right)^{5/2}\left[1+\frac{3}{2}\left(K_0^\prime -4\right) f\right],
   :label: V


.. math::
    K_{T}=(1+2f)^{5/2}\biggl[ & K_0+(3K_0{K}^\prime_{0}-5K_0)f\\ &+\frac{27}{2}(K_0{K}^\prime_{0}-4K_0)f^2 \biggr],
    :label: K


.. math::
	G = (1+& 2f)^{5/2}  \biggl[G_0+(3K_0{G}^\prime_{0}-5G_0)f\\ & +(6K_0{G}^\prime_{0}-24K_0-14G_{0}
	 + \frac{9}{2}K_{0}{K}^\prime_{0})f^2 \biggr].
    :label: G

Here :math:`K_0` and :math:`G_0` are the reference bulk modulus and shear modulus and :math:`K_0^\prime` and :math:`{G}^\prime_{0}` are the derivative of the respective moduli with respect to pressure.

BurnMan has the option to use the second-order expansion for shear modulus by dropping the :math:`f^2` terms in these equations (as is sometimes done for experimental fits or EoS modeling).

XXX include table with parametersXXX


Thermal Corrections
^^^^^^^^^^^^^^^^^^^

Thermal corrections for  pressure, and isothermal bulk modulus and shear modulus are derived from the Mie-Grüneisen-Debye EoS with the quasi-harmonic approximation.
Here we adopt the formalism of :cite:`stixrude2005` where these corrections are added to equations :eq:`V`--:eq:`G`:

.. math::
    P_{th}(V,T) &={\frac{\gamma \Delta \mathcal{U}}{V}}, \\
    :label: Pth
    K_{th}(V,T) &=(\gamma +1-q)\frac{\gamma \Delta \mathcal{U}}{V} -\gamma ^{2} \frac{\Delta(C_{V}T)}{V} ,\\
    G_{th}(V,T) &=  -\frac{\eta_{S} \Delta \mathcal{U}}{V}.

The :math:`\Delta` refers to the difference in the relevant quantity from the reference temperature (300 K).
:math:`\gamma` is the Grüneisen parameter, :math:`q` is the logarithmic volume derivative of the Grüneisen parameter, :math:`\eta_{S}` is the shear strain derivative of the Grüneisen parameter, :math:`C_V` is the heat capacity at constant volume, and :math:`\mathcal{U}` is the internal energy at temperature :math:`T`.
 :math:`C_V` and :math:`\mathcal{U}` are calculated using the Debye model for vibrational energy of a lattice.
These quantities are calculated as follows:

.. math::
    C_V &= 9nR\left (  \frac{T}{\theta}\right )^3\int_{0}^{\frac{\theta}{T}}\frac{e^{\tau}\tau^{4}}{(e^{\tau}-1)^2}d\tau, \\
    \mathcal{U} &= 9nRT\left ( \frac{T}{\theta} \right )^3\int_{0}^{\frac{\theta}{T}}\frac{\tau^3}{(e^{\tau}-1)}d\tau, \\
    \gamma &= \frac{1}{6}\frac{\nu_{0}^2}{\nu^{2}}(2f+1)\left [  a_{ii}^{(1)} +a_{iikk}^{(2)}f\right ], \\
    q &= \frac{1}{9\gamma}\left [ 18\gamma^{2}-6\gamma -\frac{1}{2} \frac{\nu^{2}_0}{\nu^2}(2f+1)^{2}a_{iikk}^{(2)} \right ], \\
    \eta_S &=-\gamma-\frac{1}{2}\frac{\nu_{0}^2}{\nu^2}(2f+1)^{2}a_{S}^{(2)}, \\
    \frac{\nu^2}{\nu^2_0} &= 1+a_{ii}^{(1)}f+\frac{1}{2}a_{iikk}^{(2)}f^2, \\
    a_{ii}^{(1)} &= 6\gamma _0, \\
    a_{iikk}^{(2)} &= -12\gamma _0+36\gamma_{0}^{2}-18q_{0}\gamma_0,  \\
    a_{S}^{(2)} &=-2\gamma _0-2\eta_{S0},

where :math:`\theta` is the Debye temperature of the mineral, :math:`\nu` is the frequency of vibrational modes for the mineral, :math:`n` is the number of atoms per formula unit (e.g. 2 for periclase, 5 for perovskite), and :math:`R` is the gas constant.
Under the approximation that the vibrational frequencies behave the same under strain, we may identify :math:`\nu/\nu_0 = \theta/\theta_0`.
The quantities :math:`\gamma_0`, :math:`\eta_{S0}` :math:`q_0`, and :math:`\theta_0` are the experimentally determined values for those parameters at the reference state.


Due to the fact that a planetary mantle is rarely isothermal along a geotherm, It is more appropriate to use the adiabatic bulk modulus :math:`K_S` instead of :math:`K_T`, which is calculated using

.. math::
    K_S=K_{T}(1+\gamma \alpha T),
    :label: K_s

where :math:`\alpha` is the coefficient of thermal expansion


.. math::
    \alpha=\frac{\gamma C_{V}V}{K_T}.
    :label: Cv


There is no difference between the isothermal and adiabatic shear moduli for an isotropic solid.
All together this makes an eleven parameter EoS model, which is summarized in Table~ :tab:`param`.
For more details on the EoS, we refer readers to :cite:`stixrude2005`.

Calculating multi-phase seismic velocities
------------------------------------------

:label: ave
Averaging schemes
^^^^^^^^^^^^^^^^^


After the thermoelastic parameters (:math:`K_S`, :math:`G`, :math:`\rho`) of each phase are determined at each pressure and/or
temperature step, these values must be combined to determine the seismic velocity of a multiphase assemblage.
We define the volume fraction of the individual minerals in an assemblage:

.. math::
    \nu_i = n_i \frac{V_i}{V},

where :math:`V_i` and :math:`n_i` are the molar volume and the molar fractions of the :math:`i` th individual phase, and :math:`V` is the total molar volume of the assemblage:



.. math::
    V = \sum_i n_i  V_i.
    :label: composite_volume


\noindent The density of the multiphase assemblage is then


.. math::
    \rho = \sum_i \nu_i \rho_i = \frac{1}{V}\sum_i {n_i \mu_i},
    :label: composite_density

where :math:`\rho_i` is the density and :math:`\mu_i` is the molar mass of the :math:`i` th phase.


Unlike density and volume, there is no straightforward way to average the bulk and shear moduli of a multiphase rock, as it depends on the specific distribution and orientation of the constituent minerals.
BurnMan allows several schemes for averaging the elastic moduli: the Voigt and Reuss bounds, the Hashin-Shtrikman bounds, the Voigt-Reuss-Hill average, and the Hashin-Shtrikman average :cite:`watt1976`.


The Voigt average, assuming constant strain across all phases, is defined as

.. math::
    X_V = \sum_i \nu_i X_i,
    :label: voigt

where :math:`X_i` is the bulk or shear modulus for the :math:`i` th phase.
The Reuss average, assuming constant stress across all phases, is defined as

.. math::
    X_R = \left(\sum_i \frac{\nu_i}{X_i} \right)^{-1}.
    :label: reuss

The Voigt-Reuss-Hill average is the arithmetic mean of Voigt and Reuss bounds:

.. math::
    X_{VRH} = \frac{1}{2} \left( X_V + X_R \right).
    :label: vrh

The Hashin-Shtrikman bounds make an additional assumption that the distribution of the phases is statistically isotropic, and are usually much narrower than the Voigt and Reuss bounds :cite:`{watt1976}.
This may be a poor assumption in regions of Earth with high anisotropy, such as the lowermost mantle, though they are rather more physically motivated than the commonly-used Voigt-Reuss-Hill average.
In most instances, the Voigt-Reuss-Hill average and the arithmetic mean of the Hashin-Shtrikman bounds are quite close to each other with the pure arithmetic mean (linear averaging) being well outside of both Hashin-Shtrikman and Voigt-Reuss-Hill.

It is worth noting that each of the above bounding methods are derived from mechanical models of a linear elastic composite.
It is thus only appropriate to apply them to elastic moduli, and not to other thermoelastic properties, such as wave speeds or density.



Computing seismic velocities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the moduli for the multiphase assemblage are computed, the compressional (:math:`P`), shear (:math:`S`) and bulk sound (:math:`\Phi`)
velocities are then result from the equations:


.. math::
    V_P = \sqrt{ \frac{K_S + \frac{4}{3} G} {\rho} }, \qquad
    V_S = \sqrt{ \frac{G}{\rho} }, \qquad
    V_\Phi = \sqrt{ \frac{K_S}{\rho} }.
    :label: seismic

To correctly compare to observed seismic velocities one needs to correct for the frequency sensitivity of attenuation.
Moduli parameters are obtained from experiments that are done at high frequencies (MHz-GHz) compared to seismic frequencies (mHz-Hz).
The frequency sensitivity of attenuation causes slightly lower velocities for seismic waves than they would be for high frequency waves.
In BurnMan one can correct the calculated acoustic velocity values to those for long period seismic tomography :cite:`Minster1981`:

.. math::
    V_{S/P}=V_{S/P}^{\mathrm{uncorr.}}\left(1-\frac{1}{2}\cot(\frac{\beta\pi}{2})\frac{1}{Q_{S/P}}(\omega)\right).

Similar to :cite:`matas2007`, we use a :math:`\beta` value of 0.3, which falls in the range of values of :math:`0.2` to :math:`0.4` proposed for the lower mantle (e.g. :cite:`karato1990`).
The correction is implemented for Q values of PREM for the lower mantle.
As :math:`Q_S` is smaller than :math:`Q_P`, the correction is more significant for S waves.
In both cases, though, the correction is minor compared to, for example, uncertainties in the temperature (corrections) and mineral physical parameters.
More involved models of relaxation mechanisms can be implemented, but lead to the inclusion of more poorly constrained parameters, :cite:`matas2007a`.
While attenuation can be ignored in many applications :cite:`trampert2001`, it might play a significant role in explaining strong variations in seismic velocities in the lowermost mantle :cite:`davies2012`.


.. _input:
User input
----------



Mineralogical composition
^^^^^^^^^^^^^^^^^^^^^^^^^

A number of pre-defined minerals are included in the mineral library and users can create their own.
The library includes wrapper functions to include a transition from the high-spin mineral to the low-spin mineral :cite:`[review:][]{lin2013} or to combine minerals for a given iron number.


*Standard minerals* -- The 'standard' mineral format includes a list of parameters given in Table :tab:`param`.
Each mineral includes a suggested EoS with which the mineral parameters are derived.
For some minerals the parameters for the thermal corrections are not yet measured or calculated, and therefore the corrections can not be applied.
An occasional mineral will not have a measured or calculated shear moduli, and therefore can only be used to compute densities and bulk sound velocities.
The mineral library is subdivided by citation.
BurnMan includes the option to produce a \LaTeX\;  table of the mineral parameters used.
BurnMan can be easily setup to incorporate uncertainties for these parameters.
*Minerals with a spin transition* -- A standard mineral for the high spin and low spin must be defined separately.
These minerals are "wrapped," so as to switch from the high spin to high spin mineral at a give pressure.
While not realistic, for the sake of simplicity, the spin transitions are considered to be sharp at a given pressure.

*Minerals depending on Fe partitioning* -- The wrapper function can partition iron, for example between ferropericlase, fp, and perovskite, pv.
It requires the input of the iron mol fraction with regards to Mg, :math:`X_\mathrm{fp}` and :math:`X_\mathrm{pv}`, which then defines the chemistry of an Mg-Fe solid solution according to (:math:`\mathrm{Mg}_{1-X_{\mathrm{Fe}}^{\mathrm{fp}}}$,$\mathrm{Fe}_{X_{\mathrm{Fe}}^{\mathrm{fp}}}$)$\mathrm{O}$ or ($\mathrm{Mg}_{1-X_{\mathrm{Fe}}^{\mathrm{pv}}}$,$\mathrm{Fe}_{X_{\mathrm{Fe}}^{\mathrm{pv}}}$)$\mathrm{SiO_3}`.
The iron mol fractions can be set to be constant or varying with P and T as needed.
Alternatively one can calculate the iron mol fraction from the distribution coefficient :math:`K_D` defined as

.. math::
    K_{D} = \frac{X_{\mathrm{Fe}}^{\mathrm{pv}}/X_{\mathrm{Mg}}^{\mathrm{pv}}}{X_{\mathrm{Fe}}^{\mathrm{fp}}/X_{\mathrm{Mg}}^{\mathrm{fp}}}.
    :label: KD


We adopt the formalism of :cite:`nakajima2012` choosing a reference distribution coefficient :math:`K_{D0}` and standard state volume change (:math:`\Delta \upsilon^{0}`) for the Fe-Mg exchange between perovskite and ferropericlase

.. math::
    K_{D}={K_D}_0 \:\exp\left(\frac{(P_0-P)\Delta \upsilon^{0}}{RT}\right),
    :label: KD2

where :math:`R` is the gas constant and :math:`P_0` the reference pressure.
As a default, we adopt the average :math:`\Delta \upsilon^{0}` of :cite:`{nakajima2012} of :math:`2\cdot10^{-7}` :math:`^3` mol:math:`^{-1}` and suggest using their :math:`{K_D}_0` value of :math:`0.5`.


The multiphase mixture of these minerals can be built by the user in three ways: 
\begin{enumerate}
\item Molar fractions of an arbitrary number of pre-defined minerals,  for example mixing standard minerals mg\_perovskite (:math:`\mathrm{MgSiO_3}`), fe\_perovskite
(:math:`\mathrm{FeSiO_3}`), periclase (:math:`\mathrm{MgO}`) and wüstite (:math:`\mathrm{FeO}`).

\item A two-phase mixture with constant or (:math:`P,T`) varying Fe partitioning using the minerals that include Fe-dependency, for example mixing :math:`\mathrm{(Mg,Fe)SiO_3}` and :math:`\mathrm{(Mg,Fe)O}` with a pre-defined distribution coefficient.

\item Weight percents (wt\%) of (Mg, Si, Fe) and distribution coefficient (includes (P,T)-dependent Fe partitioning).
This calculation assumes
that each element is completely oxidized into its corresponding oxide mineral
(:math:`\mathrm{MgO}`, :math:`\mathrm{FeO}`, :math:`\mathrm{SiO_2}`) and then combined to form iron-bearing perovskite and
ferropericlase taking into account some Fe partition coefficient.

\end{enumerate}


.. label: geothermal
Geotherm
^^^^^^^^

Unlike the pressure, the temperature of the lower mantle is relatively unconstrained.
As elsewhere, BurnMan provides a number of built-in geotherms, as well as the ability to use user-defined temperature-depth relationships.
A geotherm in BurnMan is an object that returns temperature as a function of pressure.
Alternatively, the user could ignore the geothermal and compute elastic velocities for a range of temperatures at any give pressure.

Currently, we include geotherms published by :cite:`brown1981` and :cite:`anderson1982earth`.
Alternatively one can use an adiabatic gradient defined by the thermoelastic properties of a given mineralogical model.
For a homogeneous material, the adiabatic temperature profile is given by integrating the ordinary differential equation (ODE)

.. math::
    \left(\frac{\text{d}T}{\text{d}P}\right)_S = \frac{\gamma T}{K_S}.
    :label: geoth

This equation can be extended to multiphase composite using the first law of thermodynamics to arrive at

.. math::
    \left(\frac{\text{d}T}{\text{d}P}\right)_S = \frac{ T \displaystyle\sum_{i} \frac{ n_i C_{Pi} \gamma_i }{K_{Si}}}{ \displaystyle\sum_{i} n_i C_{Pi} },
    :label: geoth2

where the subscripts correspond to the :math:`i` th phase, :math:`C_P` is the heat capacity at constant pressure of a phase, and the other symbols are as defined above.
Integrating this ODE requires a choice in anchor temperature (:math:`T_0`) at the top of the lower mantle (or including this as a parameter in an inversion).
As the adiabatic geotherm is dependent on the thermoelastic parameters at high pressures and temperatures, it is dependent on the equation of state used.


.. label: seis
Seismic Models
^^^^^^^^^^^^^^^^^^^^^^^^^


BurnMan allows for direct visual and quantitative comparison with seismic velocity models.
Various ways of plotting can be found in the examples.
Quantitative misfits between two profiles include an L2-norm and a chi-squared misfit, but user defined norms can be implemented.
A seismic model in BurnMan is
an object that provides pressure, density, and seismic velocities (:math:`V_P, V_\Phi, V_S`) as a function of depth.

To compare to seismically constrained profiles, BurnMan provides the 1D seismic velocity model PREM :cite:`{dziewonski1981}.
One can choose to evaluate :math:`V_P, V_\Phi, V_S, \rho, K_S` and/or :math:`G`.
The user can input their own seismic profile, an example of which is included for AK135 :cite:`kennett1995`.

Besides standardized 1D radial profiles, one can also compare to regionalized average profiles for the lower mantle.
This option accommodates the observation that the lowermost mantle can be clustered into two regions, a `slow' region, which represents the so-called Large Low Shear Velocity Provinces, and `fast' region, the continuous surrounding region where slabs might subduct :cite:`lekic2012`.
This clustering as well as the averaging of the 1D model occurs over five tomographic S wave velocity  models (SAW24B16: :cite:`megnin2000`; HMSL-S: :cite:`houser2008`; S362ANI: :cite:`kustowski2008`; GyPSuM: :cite:`simmons2010`; S40RTS: :cite:`ritsema2011`).
The strongest deviations from PREM occur in the lowermost 1000 km.
Using the `fast' and `slow' S wave velocity profiles is therefore most important when interpreting the lowermost mantle. Suggestion of compositional variation between these regions comes from seismology :cite:`to2005,he2012` as well as geochemistry :cite:`deschamps2012,jackson2010`.
Based on thermo-chemical convection models, :cite:`styles2011` also show that averaging profiles in thermal boundary layers may cause problems for seismic interpretation.

We additionally apply cluster analysis to and provide models for P wave velocity based on two tomographic models (MIT-P08: :cite:`li2008`; GyPSuM: :cite:`simmons2012`).
The clustering results correlate well with the fast and slow regions for S wave velocities; this could well be due to the fact that the initial model for the P wave velocity models is scaled from S wave tomographic velocity models.
Additionally, the variations in P wave velocities are a lot smaller than for S waves.
For this reason using these adapted models is most important when comparing the S wave velocities.

While interpreting lateral variations of seismic velocity in terms of composition and temperature is a major goal :cite:`trampert2004,mosca2012`, to determine the bulk composition the current challenge appears to be concurrently fitting absolute P and S wave velocities and incorporate the significant uncertainties in mineral physical parameters).


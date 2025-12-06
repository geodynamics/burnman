.. _ref-materials-mineral-equations-of-state:

Equations of State for Minerals
===============================

To calculate the bulk (:math:`K`) modulus, shear modulus (:math:`G`) and
density (:math:`\rho`) of a material at a given pressure (:math:`P`) and
temperature (:math:`T`), optionally defined by a geotherm) and determine the
seismic velocities (:math:`V_S, V_P, V_\Phi`), one uses an Equation of State
(EoS).  Currently the following EoSs are supported in BurnMan:

* Birch-Murnaghan finite-strain EoS (excludes temperature effects, :cite:`Poirier1991`), 
* Birch-Murnaghan finite-strain EoS with a Mie-Grüneisen-Debye thermal correction, as formulated by :cite:`Stixrude2005`.
* Birch-Murnaghan finite-strain EoS with a Mie-Grüneisen-Debye thermal correction, as formulated by :cite:`Matas2007`.
* Modified Tait EoS (excludes temperature effects, :cite:`HC1974`), 
* Modified Tait EoS with a pseudo-Einstein model for thermal corrections, as formulated by :cite:`HP2011`.
* Compensated-Redlich-Kwong for fluids, as formulated by :cite:`HP1991`.


To calculate these thermoelastic parameters, the EoS
requires the user to input the pressure, temperature, and the phases
and their molar fractions.  These inputs and outputs are further discussed in
:ref:`ref-methods-user-input`.



Birch-Murnaghan (isothermal)
""""""""""""""""""""""""""""

The Birch-Murnaghan equation is an isothermal Eulerian finite-strain EoS
relating pressure and volume.  The negative finite-strain (or compression) is
defined as

.. math::
    f=\frac{1}{2}\left[\left(\frac{V}{V_0}\right)^{-2/3}-1\right],
    :label: f

where :math:`V` is the volume at a given pressure and :math:`V_0` is the
volume at a reference state (:math:`P = 10^5` Pa , :math:`T` = 300 K).  The
pressure and elastic moduli are derived from a third-order Taylor expansion of
Helmholtz free energy in :math:`f` and evaluating the appropriate volume and
strain derivatives (e.g., :cite:`Poirier1991`).  For an isotropic
material one obtains for the pressure, isothermal bulk modulus, and shear
modulus:


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

Here :math:`K_0` and :math:`G_0` are the reference bulk modulus and shear
modulus and :math:`K_0^\prime` and :math:`{G}^\prime_{0}` are the derivative
of the respective moduli with respect to pressure.

BurnMan has the option to use the second-order expansion for shear modulus by
dropping the :math:`f^2` terms in these equations (as is sometimes done for
experimental fits or EoS modeling).

Modified Tait (isothermal)
""""""""""""""""""""""""""

The Modified Tait equation of state was developed by :cite:`HC1974`. It has the considerable benefit of allowing volume to be expressed as a function of pressure. It performs very well to pressures and temperatures relevant to the deep Earth :cite:`HP2011`.

.. math::
    \frac{V_{P, T}}{V_{1 bar, 298 K}} &= 1 - a(1-(1+bP)^{-c}), \\
    a &= \frac{1 + K_0'}{1 + K_0' + K_0K_0''}, \\
    b &= \frac{K_0'}{K_0} - \frac{K_0''}{1 + K_0'}, \\
    c &= \frac{1 + K_0' + K_0K_0''}{K_0'^2 + K_0' - K_0K_0''}
    :label: mtait


Mie-Grüneisen-Debye (thermal correction to Birch-Murnaghan)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Debye model for the Helmholtz free energy can be written as follows :cite:`Matas2007`

.. math::
    \mathcal{F} &= \frac{9nRT}{V}\frac{1}{x^3} \int_{0}^{x} \xi^2 \ln (1-e^{-\xi}) d\xi, \\
    x &= \theta / T, \\
    \theta &= \theta_0 \exp \left( \frac{\gamma_0-\gamma}{q_0} \right), \\
    \gamma &= \gamma_0 \left( \frac{V}{V_0} \right)^{q_0}

where :math:`\theta` is the Debye temperature and :math:`\gamma` is the Grüneisen parameter. 

Using thermodynamic relations we can derive equations for the thermal pressure and bulk modulus

.. math::
    P_{th}(V,T) &= - \frac{\partial \mathcal{F(V, T)}}{\partial V}, \\
    &= \frac{3 n \gamma R T}{V} D(x), \\
    K_{th}(V,T) &= -V \frac{\partial P(V, T)}{\partial V}, \\
    &= \frac{3 n \gamma R T}{V} \gamma \left[ (1-q_0 - 3 \gamma) D(x) + 3\gamma \frac{x}{e^x - 1} \right], \\
    D(x) &= \frac{3}{x^3} \int_{0}^{x} \frac{\xi^3}{e^{\xi} - 1} d\xi

The thermal shear correction used in BurnMan was developed by :cite:`HS1998`

.. math::
    G_{th}(V,T) = \frac{3}{5} \left[ K_{th} (V, T) - 2\frac{3nRT}{V}\gamma D(x) \right]


The total pressure, bulk and shear moduli can be calculated from the following sums

.. math::
    P(V, T) &= P_{\textrm{ref}}(V, T_0) + P_{th}(V, T) - P_{th}(V, T_0), \\
    K(V, T) &= K_{\textrm{ref}}(V, T_0) + K_{th}(V, T) - K_{th}(V, T_0), \\
    G(V, T) &= G_{\textrm{ref}}(V, T_0) + G_{th}(V, T) - G_{th}(V, T_0) 

This equation of state is substantially the same as that in SLB2005 (see below).  
The primary differences are in the thermal correction to the shear modulus and in the volume dependences of the Debye temperature and the Gruneisen parameter.


HP2011 (thermal correction to Modified Tait)
""""""""""""""""""""""""""""""""""""""""""""

The thermal pressure can be incorporated into the Modified Tait equation of state, replacing :math:`P` with :math:`P-\left(P_{\textrm{th}} - P_{\textrm{th0}}\right)` in Equation :eq:`mtait` :cite:`HP2011`. Thermal pressure is calculated using a Mie-Grüneisen equation of state and an Einstein model for heat capacity, even though the Einstein model is not actually used for the heat capacity when calculating the enthalpy and entropy (see following section).

.. math::
    P_{\textrm{th}} &= \frac{\alpha_0 K_0 E_{\textrm{th}} }{C_{V0}}, \\
    E_{\textrm{th}} &= 3 n R \Theta \left(0.5 + \frac{1}{ \exp(\frac{\Theta}{T}) - 1 }\right), \\
    C_{V} &= 3 n R \frac{(\frac{\Theta}{T})^2\exp(\frac{\Theta}{T})}{(\exp(\frac{\Theta}{T})-1)^2}


:math:`\Theta` is the Einstein temperature of the crystal in Kelvin, approximated for a substance :math:`i` with :math:`n_i` atoms in the unit formula and a molar entropy :math:`S_i` using the empirical formula 

.. math::
    \Theta_i=\frac{10636}{S_i/n_i + 6.44}


SLB2005 (for solids, thermal)
"""""""""""""""""""""""""""""

Thermal corrections for pressure, and isothermal bulk modulus and shear
modulus are derived from the Mie-Grüneisen-Debye EoS with the quasi-harmonic
approximation.  Here we adopt the formalism of :cite:`Stixrude2005` where
these corrections are added to equations :eq:`V`--:eq:`G`:

.. math::
    P_{th}(V,T) &={\frac{\gamma \Delta \mathcal{U}}{V}}, \\
    K_{th}(V,T) &=(\gamma +1-q)\frac{\gamma \Delta \mathcal{U}}{V} -\gamma ^{2} \frac{\Delta(C_{V}T)}{V} ,\\
    G_{th}(V,T) &=  -\frac{\eta_{S} \Delta \mathcal{U}}{V}.
    :label: Pth

The :math:`\Delta` refers to the difference in the relevant quantity from the
reference temperature (300 K).  :math:`\gamma` is the Grüneisen parameter,
:math:`q` is the logarithmic volume derivative of the Grüneisen parameter,
:math:`\eta_{S}` is the shear strain derivative of the Grüneisen parameter,
:math:`C_V` is the heat capacity at constant volume, and :math:`\mathcal{U}`
is the internal energy at temperature :math:`T`.  :math:`C_V` and
:math:`\mathcal{U}` are calculated using the Debye model for vibrational
energy of a lattice. These quantities are calculated as follows:

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

where :math:`\theta` is the Debye temperature of the mineral, :math:`\nu` is
the frequency of vibrational modes for the mineral, :math:`n` is the number of
atoms per formula unit (e.g. 2 for periclase, 5 for perovskite), and :math:`R`
is the gas constant.  Under the approximation that the vibrational frequencies
behave the same under strain, we may identify :math:`\nu/\nu_0 =
\theta/\theta_0`.  The quantities :math:`\gamma_0`, :math:`\eta_{S0}`
:math:`q_0`, and :math:`\theta_0` are the experimentally determined values for
those parameters at the reference state.


Due to the fact that a planetary mantle is rarely isothermal along a geotherm,
it is more appropriate to use the adiabatic bulk modulus :math:`K_S` instead
of :math:`K_T`, which is calculated using

.. math::
    K_S=K_{T}(1+\gamma \alpha T),
    :label: K_s

where :math:`\alpha` is the coefficient of thermal expansion:


.. math::
    \alpha=\frac{\gamma C_{V}V}{K_T}.
    :label: Cv


There is no difference between the isothermal and adiabatic shear moduli for
an isotropic solid.  All together this makes an eleven parameter EoS model,
which is summarized in the Table below. For more details on the
EoS, we refer readers to :cite:`Stixrude2005`.

.. _table-parameters:

+--------------+------------------+-----------------------------------+-------------------------+
|User Input    |Symbol            |Definition                         |Units                    |
|              |                  |                                   |                         |
+==============+==================+===================================+=========================+
|V_0           |:math:`V_{0}`     |Volume at P = :math:`10^5`         |m :math:`^{3}`           |
|              |                  | Pa , T = 300 K                    |mol :math:`^{-1}`        |
+--------------+------------------+-----------------------------------+-------------------------+
|K_0           |:math:`K_{0}`     |Isothermal bulk modulus at `P=10^5`|Pa                       |
|              |                  |Pa, T = 300 K                      |                         |
+--------------+------------------+-----------------------------------+-------------------------+
|Kprime_0      |:math:`K^\prime_0`|Pressure derivative of             |                         |
|              |                  |:math:`K_{0}`                      |                         |
|              |                  |                                   |                         |
+--------------+------------------+-----------------------------------+-------------------------+
|G_0           |:math:`G_{0}`     |Shear modulus at P = :math:`10^5`  |Pa                       |
|              |                  |Pa, T = 300 K                      |                         |
|              |                  |                                   |                         |
|              |                  |                                   |                         |
|              |                  |                                   |                         |
+--------------+------------------+-----------------------------------+-------------------------+
|Gprime_0      |:math:`G^\prime_0`|Pressure derivative of             |                         |
|              |                  |:math:`G_{0}`                      |                         |
|              |                  |                                   |                         |
+--------------+------------------+-----------------------------------+-------------------------+
|molar_mass    |:math:`\mu`       |mass per mole formula unit         |kg                       |
|              |                  |                                   |:math:`\mathrm{mol}^{-1}`|
|              |                  |                                   |                         |
|              |                  |                                   |                         |
+--------------+------------------+-----------------------------------+-------------------------+
|n             |n                 |number of atoms per formula unit   |                         |
|              |                  |                                   |                         |
|              |                  |                                   |                         |
|              |                  |                                   |                         |
+--------------+------------------+-----------------------------------+-------------------------+
|Debye_0       |:math:`\theta_{0}`|Debye Temperature                  |K                        |
|              |                  |                                   |                         |
+--------------+------------------+-----------------------------------+-------------------------+
|grueneisen_0  |:math:`\gamma_{0}`|Grüneisen parameter at P =         |                         |
|              |                  |:math:`10^5` Pa, T = 300 K         |                         |
|              |                  |                                   |                         |
|              |                  |                                   |                         |
+--------------+------------------+-----------------------------------+-------------------------+
|q0            |:math:`q_{0}`     |Logarithmic volume derivative of   |                         |
|              |                  |the Grüneisen parameter            |                         |
|              |                  |                                   |                         |
|              |                  |                                   |                         |
|              |                  |                                   |                         |
|              |                  |                                   |                         |
+--------------+------------------+-----------------------------------+-------------------------+
|eta_s_0       |:math:`\eta_{S0}` |Shear strain derivative of the     |                         |
|              |                  |Grüneisen parameter                |                         |
|              |                  |                                   |                         |
|              |                  |                                   |                         |
|              |                  |                                   |                         |
+--------------+------------------+-----------------------------------+-------------------------+

This equation of state is substantially the same as that of the Mie-Gruneisen-Debye (see above).  
The primary differences are in the thermal correction to the shear modulus and in the volume dependences of the Debye temperature and the Gruneisen parameter.

Compensated-Redlich-Kwong (for fluids, thermal)
"""""""""""""""""""""""""""""""""""""""""""""""

The CORK equation of state :cite:`HP1991` is a simple virial-type extension to the modified Redlich-Kwong (MRK) equation of state. It was designed to compensate for the tendency of the MRK equation of state to overestimate volumes at high pressures and accommodate the volume behaviour of coexisting gas and liquid phases along the saturation curve.

.. math::
    V &= \frac{RT}{P} + c_1 - \frac{c_0 R T^{0.5}}{(RT + c_1 P)(RT + 2 c_1 P)} + c_2 P^{0.5} + c_3 P, \\
    c_0 &= c_{0,0} T_c^{2.5}/P_c + c_{0,1} T_c^{1.5}/P_c T, \\
    c_1 &= c_{1,0} T_c/P_c, \\
    c_2 &= c_{2,0} T_c/P_c^{1.5} + c_{2,1}/P_c^{1.5} T, \\
    c_3 &= c_{3,0} T_c/P_c^2 + c_{3,1}/P_c^2 T



Calculating Thermodynamic Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

So far, we have concentrated on the thermoelastic properties of
minerals. There are, however, additional thermodynamic properties which are
required to describe the thermal properties such as the energy,
entropy and heat capacity.
These properties are related by the following expressions:

.. math::
    \mathcal{G}= \mathcal{E} - T\mathcal{S} + PV = \mathcal{H} - TS = \mathcal{F} + PV
    :label: gibbs

where :math:`P` is the pressure, :math:`T` is the temperature and :math:`\mathcal{E}`, :math:`\mathcal{F}`, :math:`\mathcal{H}`, :math:`\mathcal{S}` and :math:`V` are the molar internal energy, Helmholtz free energy, enthalpy, entropy and volume respectively.

HP2011
""""""

.. math::
    \mathcal{G}(P,T) &= \mathcal{H}_{\textrm{1 bar, T}} - T\mathcal{S}_{\textrm{1 bar, T}} + \int_{\textrm{1 bar}}^P V(P,T)~dP, \\
    \mathcal{H}_{\textrm{1 bar, T}} &= \Delta_f\mathcal{H}_{\textrm{1 bar, 298 K}} + \int_{298}^T C_P~dT, \\
    \mathcal{S}_{\textrm{1 bar, T}} &= \mathcal{S}_{\textrm{1 bar, 298 K}} + \int_{298}^T \frac{C_P}{T}~dT, \\
    \int_{\textrm{1 bar}}^P V(P,T)~dP &= P V_0 \left( 1 - a + \left( a \frac{(1-b P_{th})^{1-c} - (1 + b(P-P_{th}))^{1-c}}{b (c-1) P} \right) \right)
    :label: gibbs_hp2011


The heat capacity at one bar is given by an empirical polynomial fit to experimental data

.. math::
    C_p = a + bT + cT^{-2} + dT^{-0.5}

The entropy at high pressure and temperature can be calculated by differentiating the expression for :math:`\mathcal{G}` with respect to temperature

.. math::
    \mathcal{S}(P,T) &= S_{\textrm{1 bar, T}} + \frac{\partial  \int V dP }{\partial T}, \\
    \frac{\partial  \int V dP }{\partial T} &= V_0 \alpha_0 K_0 a \frac{C_{V0}(T)}{C_{V0}(T_\textrm{{ref}})} ((1+b(P-P_{th}))^{-c} - (1-bP_{th})^{-c} )

Finally, the enthalpy at high pressure and temperature can be calculated

.. math::
    \mathcal{H}(P,T) = \mathcal{G}(P,T) + T\mathcal{S}(P,T)

SLB2005
"""""""

The Debye model yields the Helmholtz free energy and entropy due to lattice vibrations

.. math::
    \mathcal{G} &= \mathcal{F} + PV, \\
    \mathcal{F} &= nRT \left(3 \ln( 1 - e^{-\frac{\theta}{T}}) - \int_{0}^{\frac{\theta}{T}}\frac{\tau^3}{(e^{\tau}-1)}d\tau \right), \\
    \mathcal{S} &= nR \left(4 \int_{0}^{\frac{\theta}{T}}\frac{\tau^3}{(e^{\tau}-1)}d\tau - 3 \ln(1-e^{-\frac{\theta}{T}}) \right), \\
    \mathcal{H} &= \mathcal{G} + T\mathcal{S}


Property modifiers
^^^^^^^^^^^^^^^^^^

The thermodynamic models above consider the effects of strain and quasiharmonic
lattice vibrations on the free energies of minerals at given temperatures and
pressures. There are a number of additional processes, such as isochemical
order-disorder and magnetic effects which also contribute to the total free energy
of a phase. Corrections for these additional processes can be applied in a number
of different ways. Burnman currently includes implementations of the following:

* Linear excesses (useful for DQF modifications for :cite:`HP2011`)
* Tricritical Landau model (two formulations)
* Bragg-Williams model
* Magnetic excesses

In all cases, the excess Gibbs free energy :math:`\mathcal{G}` and first and second partial derivatives with
respect to pressure and temperature are calculated. The thermodynamic properties of each
phase are then modified in a consistent manner; specifically:

.. math::
    \mathcal{G} = \mathcal{G}_o + \mathcal{G}_m, \\
    \mathcal{S} = \mathcal{S}_o - \frac{\partial \mathcal{G}}{\partial T}_m, \\
    \mathcal{V} = \mathcal{V}_o + \frac{\partial \mathcal{G}}{\partial P}_m, \\
    K_T = \mathcal{V} / \left( \frac{\mathcal{V}_o}{K_To} - \frac{\partial^2\mathcal{G}}{\partial P^2} \right)_m, \\
    C_p = C_po - T \frac{\partial ^2 \mathcal{G}}{\partial T^2}_m, \\
    \alpha = \left( \alpha_o \mathcal{V}_o  + \frac{\partial^2 \mathcal{G}}{\partial P \partial T}_m \right) / \mathcal{V}, \\
    \mathcal{H} = \mathcal{G} + T \mathcal{S}, \\
    \mathcal{F} = \mathcal{G} - P \mathcal{V}, \\
    C_v = C_p - \mathcal{V} T \alpha^2 K_T, \\
    \gamma = \frac{\alpha K_T \mathcal{V}}{C_v}, \\
    K_S = K_T \frac{C_p}{C_v}

Subscripts :math:`_o` and :math:`_m` indicate original properties and modifiers respectively.
Importantly, this allows us to stack modifications such as multiple Landau transitions
in a simple and straightforward manner. In the burnman code, we add property modifiers as an attribute to each
mineral as a list. For example::

    from burnman.minerals import SLB_2011
    stv = SLB_2011.stishovite()
    stv.property_modifiers = [
    ['landau',
    {'Tc_0': -4250.0, 'S_D': 0.012, 'V_D': 1e-09}]]
    ['linear',
    {'delta_E': 1.e3, 'delta_S': 0., 'delta_V': 0.}]]

Each modifier is a list with two elements, first the name of the modifier type, and second a dictionary with
the required parameters for that model. A list of parameters for each model is given in the following sections.

Linear excesses (linear)
""""""""""""""""""""""""
A simple linear correction in pressure and temperature. Parameters are 'delta_E', 'delta_S' and 'delta_V'.

.. math::
   \mathcal{G} = \Delta \mathcal{E} - T \Delta \mathcal{S} + P \Delta \mathcal{V}, \\
   \frac{\partial \mathcal{G}}{\partial T} = - \Delta \mathcal{S}, \\
   \frac{\partial \mathcal{G}}{\partial P} = \Delta \mathcal{V}, \\
   \frac{\partial^2 \mathcal{G}}{\partial T^2} = 0, \\
   \frac{\partial^2 \mathcal{G}}{\partial P^2} = 0, \\
   \frac{\partial^2 \mathcal{G}}{\partial T \partial P} = 0

Tricritical Landau model (landau)
"""""""""""""""""""""""""""""""""

Applies a tricritical Landau correction to the properties
of an endmember which undergoes a displacive phase transition. These
transitions are not associated with an activation energy, and therefore
they occur rapidly compared with seismic wave propagation. Parameters are
'Tc_0', 'S_D' and 'V_D'.

This correction follows :cite:`Putnis1992`, and is done relative to
the completely *ordered* state (at 0 K).
It therefore differs in implementation from both :cite:`Stixrude2011` and
:cite:`HP2011`, who compute properties relative to
the completely disordered state and standard states respectively.
The current implementation is preferred, as the excess
entropy (and heat capacity) terms are equal to zero at 0 K.

.. math::
   Tc = Tc_0 + \frac{V_D P}{S_D}

If the temperature is above the critical temperature,
Q (the order parameter) is equal to zero, and the Gibbs free energy is simply that of the disordered phase:

.. math::
   \mathcal{G}_{\textrm{dis}} = -S_D \left( \left( T - Tc \right) + \frac{Tc_0}{3} \right), \\
   \frac{\partial \mathcal{G}}{\partial P}_{\textrm{dis}} = V_D, \\
   \frac{\partial \mathcal{G}}{\partial T}_{\textrm{dis}} = -S_D

If temperature is below the critical temperature, Q is between 0 and 1.
The gibbs free energy can be described thus:

.. math::
   Q^2 = \sqrt{\left( 1 - \frac{T}{Tc} \right)}, \\
   \mathcal{G} = S_D \left((T - Tc) Q^2 + \frac{Tc_0 Q^6}{3} \right) + \mathcal{G}_{\textrm{dis}}, \\
   \frac{\partial \mathcal{G}}{\partial P} = - V_D Q^2 \left(1 + \frac{T}{2 Tc} \left(1. - \frac{Tc_0}{Tc} \right) \right) + \frac{\partial \mathcal{G}}{\partial P}_{\textrm{dis}}, \\
   \frac{\partial \mathcal{G}}{\partial T} = S_D Q^2 \left(\frac{3}{2} - \frac{Tc_0}{2 Tc} \right) + \frac{\partial \mathcal{G}}{\partial T}_{\textrm{dis}}, \\
   \frac{\partial^2 \mathcal{G}}{\partial P^2} = V_D^2 \frac{T}{S_D Tc^2 Q^2} \left( \frac{T}{4 Tc} \left(1. + \frac{Tc_0}{Tc} \right) + Q^4 \left(1. - \frac{Tc_0}{Tc} \right) - 1 \right), \\
   \frac{\partial^2 \mathcal{G}}{\partial T^2} = - \frac{S_D}{Tc Q^2} \left(\frac{3}{4} - \frac{Tc_0}{4 Tc} \right), \\
   \frac{\partial^2 \mathcal{G}}{\partial P \partial T} = \frac{V_D}{2 Tc Q^2} \left(1 + \left(\frac{T}{2 Tc} - Q^4 \right) \left(1 - \frac{Tc_0}{Tc} \right) \right)



Tricritical Landau model (landau_hp)
""""""""""""""""""""""""""""""""""""

Applies a tricritical Landau correction similar to that described
above. However, this implementation follows :cite:`HP2011`,
who compute properties relative to the standard state. Parameters are
'P_0', 'T_0', 'Tc_0', 'S_D' and 'V_D'.

It is worth noting that the correction described by :cite:`HP2011`
has been incorrectly used throughout the geological literature,
particularly in studies involving magnetite (which includes studies
comparing oxygen fugacities to the FMQ buffer (due to an incorrect
calculation of the properties of magnetite). Note that even if the
implementation is correct, it still allows the order parameter Q to be
greater than one, which is physically impossible.

We include this implementation in order to reproduce the dataset of
:cite:`HP2011`. If you are creating your own minerals, we recommend
using the standard implementation.

.. math::
   Tc = Tc0 + \frac{V_D P}{S_D}

If the temperature is above the critical temperature,
Q (the order parameter) is equal to zero. Otherwise

.. math::
   Q^2 = \sqrt{\left( \frac{Tc - T}{Tc0} \right)}

.. math::
   \mathcal{G} = Tc_0 S_D \left( Q_0^2 - \frac{Q_0 ^ 6}{3} \right) \
        - S_D \left( Tc Q^2 - Tc_0 \frac{Q ^ 6}{3} \right) \
        - T S_D \left( Q_0^2 - Q^2 \right) + P V_D Q_0^2, \\
   \frac{\partial \mathcal{G}}{\partial P} = -V_D \left( Q^2 - Q_0^2
   \right), \\
   \frac{\partial \mathcal{G}}{\partial T} = S_D \left( Q^2 - Q_0^2
   \right), \\

The second derivatives of the Gibbs free energy are only non-zero if
the order parameter exceeds zero. Then

.. math::
        \frac{\partial^2 \mathcal{G}}{\partial P^2}  = -\frac{V_D^2}{2
        S_D Tc_0 Q^2}, \\
        \frac{\partial^2 \mathcal{G}}{\partial T^2}  =
        -\frac{S_D}{2 Tc_0 Q^2}, \\
        \frac{\partial^2 \mathcal{G}}{\partial P \partial T} =
        \frac{V_D}{2 Tc_0 Q^2}


Bragg-Williams model (bragg_williams)
"""""""""""""""""""""""""""""""""""""
The Bragg-Williams model is a symmetric solution
model between endmembers with an excess configurational entropy term
determined by the specifics of order-disorder in the
mineral, multiplied by some empirical factor. Expressions for the
excess Gibbs free energy can be found in :cite:`HP1996`. Parameters are
'deltaH', 'deltaV', 'Wh', 'Wv', 'n' and 'factor'.

Magnetic model (magnetic_chs)
"""""""""""""""""""""""""""""

This model approximates the excess energy due to magnetic ordering. It
was originally described in :cite:`CHS1987`. The expressions used
by BurnMan can be found in :cite:`Sundman1991`. Parameters are
'structural_parameter', 'curie_temperature'[2] (zero pressure value
and pressure dependence) and 'magnetic_moment'[2] (zero pressure value
and pressure dependence).

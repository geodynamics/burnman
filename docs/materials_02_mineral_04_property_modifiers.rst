.. _ref-property-modifiers:

Mineral Property Modifiers
~~~~~~~~~~~~~~~~~~~~~~~~~~

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

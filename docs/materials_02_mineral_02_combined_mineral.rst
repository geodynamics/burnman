.. _ref-combined-mineral:

Combined Mineral Class
~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: burnman

The :class:`CombinedMineral` class is a subclass of :class:`Mineral`
that allows combining multiple mineral phases into a single effective phase.
This is useful, for example, when there is not enough information to
constrain a complete equation of state for a mineral. 

The essential idea is to represent the desired mineral as a weighted
mechanical mixture of other minerals with well-defined equations
of state. In practice, it is usually desirable to combine minerals
that are similar in structure and composition to the target mineral.
If the composite minerals are picked carefully enough, it may be
that only a linear correction to the Gibbs energy
of the combined phase (as a function of pressure and temperature)
is sufficient to reproduce experimental data for the target mineral.

The :class:`CombinedMineral` class takes a list of
component minerals and their corresponding molar amounts,
along with optional energy adjustments,
and computes the overall thermodynamic properties self-consistently.
The linear adjustments to the Gibbs energy are given in the order:
(1) internal energy [J/mol], (2) entropy [J/(mol·K)],
and (3) volume [m³/mol]:

.. code-block:: python

    from burnman import CombinedMineral

    component_minerals = [mineral1, mineral2, mineral3]
    molar_amounts = [0.5, 0.3, 0.2]
    energy_adjustments = [1000.0, -5.0, 1.e-6]

    combined_mineral = CombinedMineral(
        component_minerals,
        molar_amounts,
        energy_adjustments=energy_adjustments
    )

In practice, it is usually desirable to combine minerals that are
similar in structure and composition to the target mineral. If
the composite minerals are picked carefully enough, it may be
that only a linear correction to the Gibbs energy
of the combined phase (as a function of pressure and temperature)
is sufficient to reproduce experimental data for the target mineral.
More complex modifications to the Gibbs energy can be
achieved by applying :ref:`ref-property-modifiers` to one or more
of the component minerals.

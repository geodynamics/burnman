# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_spintransition_thermal
------------------------------

This example illustrates how to create a non-ideal solution model
for (Mg,Fe\ :sup:`HS`\ ,Fe\ :sup:`LS`\ )O ferropericlase that has a gradual
spin transition at finite temperature.
First, we define the MgO endmember and two endmembers for the low and high spin
states of FeO. Then we create a regular/symmetric solution that
incorporates all three endmembers. The modified solution class
contains a method called set_equilibrium_composition, which calculates
the equilibrium proportions of the low and high spin phases at the
desired bulk composition, pressure and temperature.

In this example, we neglect the elastic component of mixing. We
also implicitly apply the Bragg-Williams approximation (i.e., we assume
that there is no short-range order by only incorporating interactions that are
a function of the average occupancy of species on each distinct site).
Furthermore, the one site model [Mg,Fe\ :sup:`HS`\ ,Fe\ :sup:`LS`\ ]O
explicitly precludes long range order.

*Specifically uses:*

* :func:`burnman.Mineral`
* :func:`burnman.Solution`

*Demonstrates:*

* implementation of gradual spin transition in (Mg,Fe)O at a user-defined
  pressure and temperature
"""

from __future__ import absolute_import
from __future__ import print_function

from scipy.optimize import brentq
import numpy as np
import matplotlib.pyplot as plt


import burnman
from burnman import Mineral, minerals
from burnman.tools.chemistry import formula_mass


if __name__ == "__main__":

    """
    First, we create three Mineral objects: one for periclase,
    one for high spin wuestite and another for low spin wuestite.

    Periclase and high spin wuestite are taken directly from
    Stixrude and Lithgow-Bertelloni (2011).
    """

    periclase = burnman.minerals.SLB_2011.periclase()
    high_spin_wuestite = burnman.minerals.SLB_2011.wuestite()

    """
    In this example, we don't go to any great lengths to fit the properties
    of low spin wuestite. We use high spin wuestite as a starting point,
    subtract the spin entropy, and make the volume smaller, to agree
    approximately with the experimental data reported by Solomatova et al. (2016).
    """
    low_spin_wuestite = burnman.minerals.SLB_2011.wuestite()
    low_spin_wuestite.property_modifiers.append(['linear',
                                                 {'delta_E': 35000.,
                                                  'delta_S': -burnman.constants.gas_constant*np.log(5),
                                                  'delta_V': -1.09e-6}])

    """
    Now, we create a class derived from Solution that contains an
    additional method called set_equilibrium_composition. This method
    finds the proportions of high and low spin wuestite that satisfy
    the provided bulk composition constraint and also ensure that the
    chemical potentials of high and low spin wuestite are equal
    (i.e. that the phase is in internal equilibrium).
    """
    class ferropericlase(burnman.Solution):
        """
        Solution class for ferropericlase
        that includes a new method "set_equilibrium_composition"
        that finds the equilibrium distribution of high spin and low spin iron
        at the current state.
        """

        def __init__(self, molar_fractions=None):
            self.name = 'ferropericlase'
            self.solution_type = 'symmetric'
            self.endmembers = [[periclase, '[Mg]O'],
                               [high_spin_wuestite, '[Fehs]O'],
                               [low_spin_wuestite, '[Fels]O']]
            self.energy_interaction = [[11.e3, 11.e3],
                                       [11.e3]]
            burnman.Solution.__init__(self,
                                      molar_fractions=molar_fractions)

        def set_equilibrium_composition(self, molar_fraction_FeO):
            """
            This method finds the equilibrium fractions of
            the three endmembers of the solution for a given bulk composition
            at the current state (pressure and temperature) of the solution.
            """
            def delta_mu(p_LS):
                """
                This function calculates the difference between the
                partial gibbs energies of high and low spin FeO for a
                given proportion of low spin FeO.
                """
                self.set_composition([1. - molar_fraction_FeO,
                                      molar_fraction_FeO*(1. - p_LS),
                                      molar_fraction_FeO*p_LS])
                return self.partial_gibbs[1] - self.partial_gibbs[2]

            # Try to find the equilibrium proportion of low spin iron
            # brentq fails with a ValueError if the
            # equilibrium proportion is exactly 0 or 1,
            # in which case we find out which is more stable.
            try:
                p_LS = brentq(delta_mu, 0., 1.)
            except ValueError:
                self.set_composition([1. - molar_fraction_FeO,
                                      molar_fraction_FeO, 0.])
                G0 = self.gibbs
                self.set_composition([1. - molar_fraction_FeO, 0.,
                                      molar_fraction_FeO])
                G1 = self.gibbs
                p_LS = 0. if G0 < G1 else 1.

            # finally, set the equilibrium composition that we have
            # just calculated
            self.set_composition([1. - molar_fraction_FeO,
                                  molar_fraction_FeO*(1. - p_LS),
                                  molar_fraction_FeO*p_LS])

    # In this line, we create our solution object
    fper = ferropericlase()

    # Now we loop over a series of pressures at three different temperatures,
    # calculating the equilibrium composition of the solution at each.
    # We fix the bulk composition of the solution to be (Mg0.8Fe0.2)O.
    X_Fe = 0.2

    pressures = np.linspace(10.e9, 150.e9, 101)
    volumes = np.empty_like(pressures)
    volumes_HS = np.empty_like(pressures)
    volumes_LS = np.empty_like(pressures)
    p_LS = np.empty_like(pressures)

    fig = plt.figure(figsize=(8, 4))
    ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

    for T, color in [[300., 'blue'],
                     [1800., 'purple'],
                     [3300., 'red']]:
        for i, P in enumerate(pressures):

            # Calculate and store the equilibrium volume and proportion of
            # iron in the low spin state
            fper.set_state(P, T)
            fper.set_equilibrium_composition(X_Fe)

            volumes[i] = fper.V
            p_LS[i] = (fper.molar_fractions[2]
                       / (fper.molar_fractions[1]+fper.molar_fractions[2]))

            # Also calculate and store the volumes if all iron were in the
            # high or low spin state
            fper.set_composition([1.-X_Fe, X_Fe, 0.])
            volumes_HS[i] = fper.V
            fper.set_composition([1.-X_Fe, 0., X_Fe])
            volumes_LS[i] = fper.V

        # Do some plotting
        ax[0].fill_between(pressures/1.e9, volumes_HS*1.e6, volumes_LS
                           * 1.e6, alpha=0.15, color=color, label=f'{T} K, volume range')
        ax[0].plot(pressures/1.e9, volumes*1.e6, c=color,
                   linewidth=2, label=f'{T} K, equilibrium volume')

        ax[1].plot(pressures/1.e9, 1.-p_LS, c=color, label=f'{T} K')

    # Add legends and axis titles to the plot
    for i in range(2):
        ax[i].legend()
        ax[i].set_xlabel('Pressure (GPa)')

    ax[0].set_ylabel('Volume (cm$^3$/mol)')
    ax[1].set_ylabel('High spin fraction')

    ax[0].set_ylim(7, 13)

    # Tidy the plot and show it
    fig.set_tight_layout(True)
    plt.show()

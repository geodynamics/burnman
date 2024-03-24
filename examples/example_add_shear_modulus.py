# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2023 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_add_shear_modulus
-------------------------

This example shows how to create a new equation of state derived
from another in order to calculate something new. In this case, we
seek to create an equation of state that builds on the Holland and
Powell (2011) equation of state, adding calculations for the shear
modulus as proposed by Hacker and Abers (2004).

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.Mineral`


*Demonstrates:*

* How to create a new equation of state and instantiate a mineral using it.
* How to output thermoelastic properties from the new mineral object.

"""
from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

from burnman import Mineral
from burnman.minerals import HGP_2018_ds633
from burnman.eos.hp import HP_TMT

if __name__ == "__main__":
    # First, we create a new class, derived from HP_TMT
    # (the modified Tait equation of state derived by Holland and Powell, 2011).
    # We overwrite just a single function; the one to calculate the shear modulus.
    # Here we use the equations given by Hacker and Abers (2004).
    class HP_TMT_shear(HP_TMT):
        def shear_modulus(self, pressure, temperature, volume, params):
            x = params["V_0"] / volume
            phi = np.log(1.0 / x)
            gamma = params["Gamma_0"]
            f = 0.5 * (np.power(x, 2.0 / 3.0) - 1.0)
            f2 = np.power(1.0 + 2.0 * f, 5.0 / 2.0)
            G_T = params["shear_modulus_0"] * np.exp(-gamma * phi)
            K_PT = self.isothermal_bulk_modulus_reuss(
                pressure, temperature, volume, params
            )
            K_T = K_PT / ((1.0 - f * (5.0 - 3.0 * params["Kprime_0"])) * f2)
            G = G_T * f2 * (1.0 - f * (5.0 - 3.0 * params["Gprime_0"] * K_T / G_T))
            return G

    # We want to create a forsterite object
    # that uses this new equation of state. For this, we need values
    # for the new parameters required by that equation of state.
    # Here we create a dictionary with the values from Hacker and Abers (2004).
    shear_params = {
        "fo": {"shear_modulus_0": 81.6e9, "Gamma_0": 5.19, "Gprime_0": 1.82}
    }

    # One way to use the new equation of state is to first create a Mineral object
    # from the existing forsterite class in HGP_2018_ds633, update the parameters
    # and then use the set_method() method:
    fo = HGP_2018_ds633.fo()
    fo.params.update(shear_params["fo"])
    fo.set_method(HP_TMT_shear)

    # Now we can use the modified object:
    P = 1.0e9
    T = 500.0
    fo.set_state(P, T)

    # Print some values:
    print(f"P: {P/1.e9:.2f} GPa, T: {T:.2f} K")

    print("\nProperties from modified object:")
    print(f"    Shear modulus: {fo.G/1.e9:.2f} GPa")
    print(f"    V_p: {fo.v_p/1.e3:.2f} km/s, V_s: {fo.v_s/1.e3:.2f} km/s")

    # This method works well if we only want one instantiation, but what if
    # we want to create several? One way to do this is to define a new class
    # and put the required steps into the initialisation function:

    class fo_w_shear(Mineral):
        def __init__(self):
            tmp = HGP_2018_ds633.fo()
            self.params = deepcopy(tmp.params)
            name = self.params["name"]
            self.params.update(shear_params[name])
            Mineral.__init__(self)
            self.set_method(HP_TMT_shear)

    # Now we can create new object from the modified class
    fo = fo_w_shear()
    fo.set_state(P, T)

    print("\nProperties from object created from modified class:")
    print(f"    Shear modulus: {fo.G/1.e9:.2f} GPa")
    print(f"    V_p: {fo.v_p/1.e3:.2f} km/s, V_s: {fo.v_s/1.e3:.2f} km/s")

    # We can now do anything we like with this new mineral,
    # like use the evaluate method:
    pressures = np.linspace(1.0e5, 10.0e9, 101)
    T = 1500.0
    temperatures = np.ones_like(pressures) * T

    v_p, v_s = fo.evaluate(["v_p", "v_s"], pressures, temperatures)
    plt.plot(pressures / 1.0e9, v_p / 1.0e3, label=f"$V_P$ for {fo.name} at {T} K")
    plt.plot(pressures / 1.0e9, v_s / 1.0e3, label=f"$V_S$ for {fo.name} at {T} K")
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Velocities (km/s)")
    plt.legend()
    plt.show()

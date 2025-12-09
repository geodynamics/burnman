# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_geotherms
-----------------

This example shows each of the geotherms currently possible with BurnMan.
These are:

1. Brown and Shankland, 1981 :cite:`Brown1981`
2. Anderson, 1982 :cite:`anderson1982earth`
3. Stacey, 1977 continental :cite:`stacey1977`
4. Stacey, 1977 oceanic :cite:`stacey1977`
5. Custom tabulated geotherm
6. Adiabatic geotherm for a given mineral assemblage

*Uses:*

* :class:`burnman.geotherm.BrownShankland`
* :class:`burnman.geotherm.Anderson`
* :class:`burnman.geotherm.StaceyContinental`
* :class:`burnman.geotherm.StaceyOceanic`
* :class:`burnman.geotherm.Geotherm`
* :class:`burnman.geotherm.AdiabaticGeotherm`
* :class:`burnman.Composite` for adiabat

*Demonstrates:*

* the available geotherms

"""

import numpy as np
import matplotlib.pyplot as plt


import burnman
from burnman import minerals


if __name__ == "__main__":
    # We want to evaluate several geotherms at particular pressures
    # The geotherms are defined as a function of depth, so we convert pressure
    # to depth using a seismic model, in this case PREM
    pressures = np.linspace(9.0e9, 128e9, 300)
    seismic_model = burnman.seismic.PREM()
    depths = seismic_model.depth(pressures)

    # Load the Earth-based geotherms packaged with BurnMan
    geotherm_BrownShankland = burnman.geotherm.BrownShankland()
    geotherm_Anderson = burnman.geotherm.Anderson()
    geotherm_StaceyContinental = burnman.geotherm.StaceyContinental()
    geotherm_StaceyOceanic = burnman.geotherm.StaceyOceanic()
    geotherm_Katsura = burnman.geotherm.Katsura2022()
    geotherm_Anzellini = burnman.geotherm.Anzellini2013(seismic_model)

    # Evaluate the temperatures at the desired depths
    temperature_BrownShankland = geotherm_BrownShankland.temperatures(depths)
    temperature_Anderson = geotherm_Anderson.temperatures(depths)
    temperature_StaceyContinental = geotherm_StaceyContinental.temperatures(depths)
    temperature_StaceyOceanic = geotherm_StaceyOceanic.temperatures(depths)
    temperature_Katsura = geotherm_Katsura.temperatures(depths)
    temperature_Anzellini = geotherm_Anzellini.temperatures(depths)

    # We can also define our own geotherm using the Geotherm class.
    depths_array = [0.0, 100.0e3, 660.0e3, 3000.0e3]
    temperatures_array = [300.0, 1500.0, 1900.0, 2500.0]
    geotherm_custom = burnman.geotherm.Geotherm(depths_array, temperatures_array)
    temperature_custom = geotherm_custom.temperatures(depths)

    # We can also calculate an adiabatic geotherm for a given material.
    # For this we need to define a material, for example a perovskite
    # and ferropericlase assemblage.
    molar_fraction_perovskite = 0.8
    fe_pv = 0.05
    fe_pc = 0.2
    pv = minerals.SLB_2011.mg_fe_perovskite()
    pc = minerals.SLB_2011.ferropericlase()
    pv.set_composition([1.0 - fe_pv, fe_pv, 0.0])
    pc.set_composition([1.0 - fe_pc, fe_pc])
    example_rock = burnman.Composite(
        [pv, pc], [molar_fraction_perovskite, 1.0 - molar_fraction_perovskite]
    )

    # Next, define an anchor temperature at the first depth
    # and calculate the adiabatic geotherm.
    # Let's choose 1600 K for the upper mantle at 9 GPa
    # (the first pressure in the list).
    T_0 = 1600.0
    geotherm_adiabatic = burnman.geotherm.AdiabaticGeotherm(
        example_rock, T_0=T_0, P_0=9.0e9, depth_to_pressure_model=seismic_model
    )
    temperature_adiabatic = geotherm_adiabatic.temperatures(depths)

    # Finally, we plot the results
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(pressures / 1e9, temperature_BrownShankland, label="Brown & Shankland")
    ax.plot(pressures / 1e9, temperature_Anderson, label="Anderson")
    ax.plot(pressures / 1e9, temperature_StaceyContinental, label="Stacey Continental")
    ax.plot(pressures / 1e9, temperature_StaceyOceanic, label="Stacey Oceanic")
    ax.plot(pressures / 1e9, temperature_Katsura, label="Katsura 2022")
    ax.plot(pressures / 1e9, temperature_Anzellini, label="Anzellini 2013")
    ax.plot(pressures / 1e9, temperature_custom, label="Custom Geotherm")
    ax.plot(
        pressures / 1e9,
        temperature_adiabatic,
        "-m",
        label="Adiabat with pv (80%) and fp(20%)",
    )

    ax.legend(loc="lower right")
    ax.set_xlim([8.5, 130])
    ax.set_xlabel("Pressure/GPa")
    ax.set_ylabel("Temperature")
    fig.savefig("output_figures/example_geotherm.png")
    plt.show()

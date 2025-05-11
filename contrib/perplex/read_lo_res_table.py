from __future__ import absolute_import
from __future__ import print_function

import os
import numpy as np
import matplotlib.pyplot as plt

from burnman.classes.perplex import PerplexMaterial

project_name = "iron_olivine_lo_res"
outfile = "iron_olivine_lo_res_table.dat"

outfile_path = os.path.join(os.getcwd(), project_name, outfile)

m = PerplexMaterial(outfile_path)

pressures = np.linspace(1.0e5, 10.0e9, 21)
temperatures = np.linspace(200.0, 2000.0, 21)

pp, tt = np.meshgrid(pressures, temperatures)
Cp, rho = m.evaluate(["C_p", "rho"], pp, tt)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns

rho_plot = ax[0].contourf(pp / 1.0e9, tt, rho, cmap="plasma", levels=21)
ax[0].set_title("Density (kg/m3)")
fig.colorbar(rho_plot, ax=ax[0])

cp_plot = ax[1].contourf(pp / 1.0e9, tt, Cp / m.molar_mass, cmap="viridis", levels=21)
ax[1].set_title("Cp (J/K/kg)")
fig.colorbar(cp_plot, ax=ax[1])

for i in range(2):
    ax[i].set_xlabel("Pressure (GPa)")
    ax[i].set_ylabel("Temperature (K)")

fig.set_layout_engine("tight")
plt.show()

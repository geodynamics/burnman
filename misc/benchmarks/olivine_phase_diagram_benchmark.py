# Benchmarks for the chemical potential functions
import burnman
from burnman.minerals import SLB_2011
from burnman import equilibrate

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Initialize the minerals we will use in this example.
ol = SLB_2011.mg_fe_olivine()
wad = SLB_2011.mg_fe_wadsleyite()
rw = SLB_2011.mg_fe_ringwoodite()

# Set the starting guess compositions for each of the solutions
ol.set_composition([0.90, 0.10])
wad.set_composition([0.90, 0.10])
rw.set_composition([0.80, 0.20])

T = 1673.15  # K

# First, we find the compositions of the three phases
# at the univariant.
composition = {"Fe": 0.2, "Mg": 1.8, "Si": 1.0, "O": 4.0}
assemblage = burnman.Composite([ol, wad, rw], [1.0, 0.0, 0.0])
equality_constraints = [
    ("T", T),
    ("phase_fraction", (ol, 0.0)),
    ("phase_fraction", (rw, 0.0)),
]
free_compositional_vectors = [{"Mg": 1.0, "Fe": -1.0}]

sol, prm = equilibrate(
    composition,
    assemblage,
    equality_constraints,
    free_compositional_vectors,
    verbose=False,
)

if not sol.success:
    raise Exception(
        "Could not find solution for the univariant using " "provided starting guesses."
    )

# We interrogate the stored copy of the assemblage for the pressure and
# the composition of all three phases
P_univariant = sol.assemblage.pressure
x_fa_univariant = sol.assemblage.phases[0].molar_fractions[1]
x_fwd_univariant = sol.assemblage.phases[1].molar_fractions[1]
x_frw_univariant = sol.assemblage.phases[2].molar_fractions[1]

print(
    f"Univariant point at T={T} K: P={P_univariant/1.0e9:.2f} GPa, "
    f"x_fa={x_fa_univariant:.3f}, x_fwd={x_fwd_univariant:.3f}, x_frw={x_frw_univariant:.3f}"
)

# Now we solve for the stable sections of the three binary loops
phase_loops = []
i = 0
for d in [
    [ol, wad, np.linspace(x_fa_univariant, 0.001, 20)],
    [ol, rw, np.linspace(x_fa_univariant, 0.999, 20)],
    [wad, rw, np.linspace(x_fwd_univariant, 0.001, 20)],
]:
    m1, m2, x_fe_m1 = d
    assemblage = burnman.Composite([m1, m2], [1.0, 0.0])

    # Reset the compositions of the two phases to have compositions
    # close to those at the univariant point
    m1.set_composition([1.0 - x_fwd_univariant, x_fwd_univariant])
    m2.set_composition([1.0 - x_fwd_univariant, x_fwd_univariant])

    # Also set the pressure and temperature
    assemblage.set_state(P_univariant, T)

    # Here our equality constraints are temperature,
    # the phase fraction of the second phase,
    # and we loop over the composition of the first phase.
    equality_constraints = [
        ("T", T),
        (
            "phase_composition",
            (m1, [["Mg_A", "Fe_A"], [0.0, 1.0], [1.0, 1.0], x_fe_m1]),
        ),
        ("phase_fraction", (m2, 0.0)),
    ]

    sols, prm = equilibrate(
        composition,
        assemblage,
        equality_constraints,
        free_compositional_vectors,
        verbose=False,
    )

    # Process the solutions
    out = np.array(
        [
            [
                sol.assemblage.pressure,
                sol.assemblage.phases[0].molar_fractions[1],
                sol.assemblage.phases[1].molar_fractions[1],
            ]
            for sol in sols
            if sol.success
        ]
    )

    phase_loops.append(out)

"""
Plot model
"""
fig1 = mpimg.imread("../../burnman/data/input_figures/slb_fig10a.png")
plt.imshow(fig1, extent=[0, 1, 0.0, 30.0], aspect="auto")

plt.plot(
    [x_fa_univariant, x_frw_univariant],
    [P_univariant / 1.0e9, P_univariant / 1.0e9],
    color="black",
    linewidth=2,
    label="univariant",
)

for loop in phase_loops:
    pressures, x_fe_m1s, x_fe_m2s = loop.T
    plt.plot(x_fe_m1s, pressures / 1.0e9, "red", linewidth=1)
    plt.plot(x_fe_m2s, pressures / 1.0e9, "red", linewidth=1)
    plt.fill_betweenx(pressures / 1.0e9, x_fe_m1s, x_fe_m2s, color="gray", alpha=0.2)

plt.title("Mg2SiO4-Fe2SiO4 phase diagram")
plt.xlabel("X_Fe")
plt.ylabel("Pressure (GPa)")
plt.legend(loc="upper right")
plt.show()

# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Pitzer_Sterner_1994
"""

import numpy as np

from burnman import Mineral
from burnman.utils.chemistry import dictionarize_formula, formula_mass


class H2O_Pitzer_Sterner(Mineral):
    def __init__(self):
        formula = "H2O"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "H2O",
            "equation_of_state": "pitzer-sterner",
            "formula": formula,
            "F_0": -249677.75163034006
            + 233.2550
            * 298.15,  # Fit to Barin at 1.e5 Pa, 2000 K, convert from HSC->SUP convention
            "Cv_0": 8.31446 * 2.25,  # Eyeball fit on Barin heat capacity
            "Debye_0": 3800,  # Eyeball fit on Barin entropy
            "Debye_n": 1.2,  # Eyeball fit on Barin entropy
            "c_coeffs": np.array(
                [
                    [0.0, 0.0, +0.24657688e6, +0.51359951e2, 0.0, 0.0],
                    [0.0, 0.0, +0.58638965e0, -0.28646939e-2, +0.31375577e-4, 0.0],
                    [
                        0.0,
                        0.0,
                        -0.62783840e1,
                        +0.14791599e-1,
                        +0.35779579e-3,
                        +0.15432925e-7,
                    ],
                    [0.0, 0.0, 0.0, -0.42719875e0, -0.16325155e-4, 0.0],
                    [0.0, 0.0, +0.56654978e4, -0.16580167e2, +0.76560762e-1, 0.0],
                    [0.0, 0.0, 0.0, +0.10917883e0, 0.0, 0.0],
                    [
                        +0.38878656e13,
                        -0.13494878e9,
                        +0.30916564e6,
                        +0.75591105e1,
                        0.0,
                        0.0,
                    ],
                    [0.0, 0.0, -0.65537898e5, +0.18810675e3, 0.0, 0.0],
                    [
                        -0.14182435e14,
                        +0.18165390e9,
                        -0.19769068e6,
                        -0.23530318e2,
                        0.0,
                        0.0,
                    ],
                    [0.0, 0.0, +0.92093375e5, +0.12246777e3, 0.0, 0.0],
                ]
            ),
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }

        Mineral.__init__(self)

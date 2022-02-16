# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Saxena and Eriksson (2015)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Iron endmember minerals and melt taken from :cite:`SE2015`
using the equation of state of :cite:`Brosh2007`.

1 bar gibbs free energy coefficients are given in the following order:
[[T_max, [const, T, T*ln(T), T^(-1), T^(-2), T^(-3), T^(-9),
T^2, T^3, T^4, T^7, T^(1/2), ln(T)]]
"""

from __future__ import absolute_import

from ..classes.mineral import Mineral
from ..utils.chemistry import formula_mass


class bcc_iron(Mineral):
    """
    BCC iron from :cite:`SE2015`.
    """

    def __init__(self):
        formula = {'Fe': 1.0}
        m = formula_mass(formula)
        self.params = {
            'name': 'BCC iron',
            'formula': formula,
            'equation_of_state': 'brosh_calphad',
            'molar_mass': m,
            'n': sum(formula.values()),
            'gibbs_coefficients': [[1811., [1225.7, 124.134, -23.5143, 77359.,
                                            0., 0., 0., -0.439752e-2,
                                            -5.8927e-8, 0., 0., 0., 0.]],
                                   [6000., [-25383.6, 299.3126, -46., 0.,
                                            0., 0., 2.29603e31, 0., 0.,
                                            0., 0., 0., 0.]],
                                   [12000., [-25383.4, 299.3122, -45.99997, 0.,
                                             0., 0., 0., 0., 0.,
                                             0., 0., 0., 0.]]],
            'V_0': 7.05e-6,  # V0
            'K_0': 170.e9,  # b8
            'Kprime_0': 6.2,  # b9
            'theta_0': 300,  # b3
            'grueneisen_0': 1.55,  # b1
            'delta': [6., 15.],  # b5, b7
            'b': [1.,  3.]  # b4, b6
        }
        self.property_modifiers = [['magnetic_chs',
                                    {'structural_parameter': 0.4,
                                     'curie_temperature': [1043., 0.],
                                     'magnetic_moment': [2.22, 0.]}]]
        Mineral.__init__(self)


class fcc_iron(Mineral):
    """
    FCC iron from :cite:`SE2015`.
    """

    def __init__(self):
        formula = {'Fe': 1.0}
        m = formula_mass(formula)
        self.params = {
            'name': 'FCC iron',
            'formula': formula,
            'equation_of_state': 'brosh_calphad',
            'molar_mass': m,
            'n': sum(formula.values()),
            'gibbs_coefficients': [[1811., [-236.7, 132.416, -24.6643, 77359.,
                                            0., 0., 0., -0.375752e-2,
                                            -5.8927e-8, 0., 0., 0., 0.]],
                                   [6000., [-27097.4, 300.2526, -46., 0.,
                                            0., 0., 2.78854e31, 0., 0.,
                                            0., 0., 0., 0.]],
                                   [12000., [-27097.1, 300.2522, -45.99996, 0.,
                                             0., 0., 0., 0., 0.,
                                             0., 0., 0., 0.]]],
            'V_0': 6.826e-6,  # V0
            'K_0': 140.e9,  # b8
            'Kprime_0': 8.,  # b9
            'theta_0': 250.,  # b3
            'grueneisen_0': 2.,  # b1
            'delta': [4., 10.],  # b5, b7
            'b': [1., 3.]  # b2, b6
        }

        self.property_modifiers = [['magnetic_chs',
                                    {'structural_parameter': 0.28,
                                     'curie_temperature': [201., 0.],
                                     'magnetic_moment': [2.1, 0.]}]]

        Mineral.__init__(self)


class hcp_iron(Mineral):
    """
    HCP iron from :cite:`SE2015`.
    """

    def __init__(self):
        formula = {'Fe': 1.0}
        m = formula_mass(formula)
        self.params = {
            'name': 'HCP iron',
            'formula': formula,
            'equation_of_state': 'brosh_calphad',
            'molar_mass': m,
            'n': sum(formula.values()),
            'gibbs_coefficients': [[1811., [-2480.08, 136.725, -24.6643,
                                            77359., 0., 0., 0., -0.375752e-2,
                                            -5.8927e-8, 0., 0., 0., 0.]],
                                   [6000., [-29340.8, 304.5616, -46., 0.,
                                            0., 0., 2.78854e31, 0., 0.,
                                            0., 0., 0., 0.]],
                                   [12000., [-29340.5, 304.5612, -45.99996, 0.,
                                             0., 0., 0., 0., 0.,
                                             0., 0., 0., 0.]]],
            'V_0': 6.677e-6,  # V0
            'K_0': 170.e9,  # b8
            'Kprime_0': 5.5,  # b9
            'theta_0': 250.,  # b3
            'grueneisen_0': 2.85,  # b1
            'delta': [6., 10.],  # b5, b7
            'b': [0.7, 2.49614]  # b4, b6
        }
        Mineral.__init__(self)


class liquid_iron(Mineral):
    """
    Liquid iron from :cite:`SE2015`.
    """

    def __init__(self):
        formula = {'Fe': 1.}
        m = formula_mass(formula)
        self.params = {
            'name': 'Liquid iron',
            'formula': formula,
            'equation_of_state': 'brosh_calphad',
            'molar_mass': m,
            'n': sum(formula.values()),
            'gibbs_coefficients': [[1811., [13265.87, 117.5756, -23.5143,
                                            77359., 0., 0., 0., -0.439752e-2,
                                            -5.8927e-8, 0., -0.3675155e-20,
                                            0., 0., 0.]],
                                   [12000., [-10838.8, 291.302, -46., 0.,
                                             0., 0., 0., 0., 0.,
                                             0., 0., 0., 0.]]],
            'V_0': 7.4602e-6,  # V0
            'K_0': 165.e9,  # b8
            'Kprime_0': 4.4729,  # b9
            'theta_0': 250.,  # b3
            'grueneisen_0': 2.,  # b1
            'delta': [6., 4.],  # b5, b7
            'b': [1., 5.10624]  # b4, b6
        }
        Mineral.__init__(self)

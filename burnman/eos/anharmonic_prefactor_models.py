# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np


class AnharmonicPrefactorModel(object):
    """
    Base class for anharmonic prefactor models.
    """

    def value(self, Vrel, params):
        raise NotImplementedError("need to implement value() in derived class!")

    def dVrel(self, Vrel, params):
        raise NotImplementedError("need to implement dVrel() in derived class!")

    def d2dVrel2(self, Vrel, params):
        raise NotImplementedError("need to implement d2dVrel2() in derived class!")


class PowerLaw(AnharmonicPrefactorModel):
    """
    Class providing methods to compute the prefactor A in the anharmonic
    contribution to thermodynamic models.

    The prefactor is defined as :math:`A = a_{anh} * (V/V_0)^{m_{anh}}`,
    with both :math:`a_{anh}` and :math:`m_{anh}` being parameters of the model.

    :return: _description_
    :rtype: _type_
    """

    def value(self, Vrel, params):
        return params["a_anh"] * np.power(Vrel, params["m_anh"])

    def dVrel(self, Vrel, params):
        return params["a_anh"] * params["m_anh"] * np.power(Vrel, params["m_anh"] - 1)

    def d2dVrel2(self, Vrel, params):
        return (
            params["a_anh"]
            * params["m_anh"]
            * (params["m_anh"] - 1)
            * np.power(Vrel, params["m_anh"] - 2)
        )


class Sigmoid(AnharmonicPrefactorModel):
    """
    Class providing methods to compute the prefactor A in the anharmonic
    contribution to thermodynamic models.

    The prefactor is defined as
    :math:`A = a_{anh} * (1 - 1/(1 + (Vrel + b_{anh}) ** c_{anh}))`
    where :math:`a_{anh}`, :math:`b_{anh}` and :math:`c_{anh}`
    are parameters of the model.

    :return: _description_
    :rtype: _type_
    """

    def value(self, Vrel, params):
        a = params["a_anh"]
        b = params["b_anh"]
        c = params["c_anh"]
        return a * (1.0 - 1.0 / (1.0 + (Vrel + b) ** c))

    def dVrel(self, Vrel, params):
        a = params["a_anh"]
        b = params["b_anh"]
        c = params["c_anh"]
        denom = (1 + (Vrel + b) ** c) ** 2
        return a * c * (Vrel + b) ** (c - 1) / denom

    def d2dVrel2(self, Vrel, params):
        a = params["a_anh"]
        b = params["b_anh"]
        c = params["c_anh"]
        num = (c - 1) - (c + 1) * (Vrel + b) ** c
        denom = (1 + (Vrel + b) ** c) ** 3
        return a * c * (Vrel + b) ** (c - 2) * num / denom

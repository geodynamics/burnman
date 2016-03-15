# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import inspect
from . import slb
from . import mie_grueneisen_debye as mgd
from . import birch_murnaghan as bm
from . import birch_murnaghan_4th as bm4
from . import modified_tait as mt
from . import hp
from . import cork
from . import vinet
from .equation_of_state import EquationOfState


def create(method):
    """
    Creates an instance of an EquationOfState from a string,
    a class EquationOfState, or an instance of EquationOfState.
    """
    if isinstance(method, str):
        if method == "slb2":
            return slb.SLB2()
        elif method == "vinet":
            return vinet.Vinet()
        elif method == "mgd2":
            return mgd.MGD2()
        elif method == "mgd3":
            return mgd.MGD3()
        elif method == "slb3":
            return slb.SLB3()
        elif method == "bm2":
            return bm.BM2()
        elif method == "bm3":
            return bm.BM3()
        elif method == "bm4":
            return bm4.BM4()
        elif method == "mt":
            return mt.MT()
        elif method == "hp_tmt":
            return hp.HP_TMT()
        elif method == "cork":
            return cork.CORK()
        else:
            raise Exception("unsupported material method " + method)
    elif isinstance(method, EquationOfState):
        return method
    elif inspect.isclass(method) and issubclass(method, EquationOfState):
        return method()
    else:
        raise Exception(
            "unsupported material method " + method.__class__.__name__)

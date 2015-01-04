# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import inspect
import slb
import mie_grueneisen_debye as mgd
import birch_murnaghan as bm
import modified_tait as mt
import hp 
import cork
from equation_of_state import EquationOfState


def create(method):
    """
    Creates an instance of an EquationOfState from a string,
    a class EquationOfState, or an instance of EquationOfState.
    """
    if isinstance(method, basestring):
        if method == "slb2":
            return slb.SLB2()
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
        elif method == "mtait":
            return mt.MT()
        elif method == "hp_tmt":
            return hp.HPMT()
        elif method == "cork":
            return cork.CORK()
        else:
            raise Exception("unsupported material method " + method)
    elif isinstance(method, EquationOfState):
        return method
    elif inspect.isclass(method) and issubclass(method, EquationOfState):
        return method()
    else:
        raise Exception("unsupported material method " + method.__class__.__name__)

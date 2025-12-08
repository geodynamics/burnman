# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


import inspect
from collections import OrderedDict
from . import slb
from . import mie_grueneisen_debye as mgd
from . import modular_mie_grueneisen_debye as mmgd
from . import murnaghan
from . import birch_murnaghan as bm
from . import modified_tait as mt
from . import macaw
from . import spock
from . import dks_liquid
from . import dks_solid
from . import hp
from . import cork
from . import vinet
from . import morse_potential
from . import reciprocal_kprime
from . import aa
from . import brosh_calphad
from .equation_of_state import EquationOfState


class CombinedMineralMethod(object):
    """Dummy class because CombinedMineral objects are derived
    from a mechanical Solution.
    Solution needs a method to call Mineral.set_state(),
    but a CombinedMineral should never have a method that
    is used for solutions."""

    def validate_parameters(self, params):
        pass

    pass


eos_names = OrderedDict()
eos_names["murnaghan"] = "Murnaghan EoS (Murnaghan, 1944)"
eos_names["bm3shear2"] = (
    "3rd-order Birch-Murnaghan with 2nd order Shear Modulus expansion (Birch, 1947)"
)
eos_names["bm3"] = "3rd-order Birch-Murnaghan EoS (Birch, 1947)"
eos_names["bm4"] = "4th-order Birch-Murnaghan EoS (Birch, 1947)"
eos_names["vinet"] = "Vinet EoS (Vinet et al., 1987)"
eos_names["morse"] = "Morse Potential EoS (Morse, 1929)"
eos_names["rkprime"] = "Reciprocal K' EoS (Stacey and Davis, 2004)"
eos_names["mgd2"] = "Mie-Gruneisen-Debye 2nd-order EoS (Matas et al., 2007)"
eos_names["mgd3"] = "Mie-Gruneisen-Debye 3rd-order EoS (Matas et al., 2007)"
eos_names["slb2"] = (
    "Stixrude and Lithgow-Bertelloni 2nd-order EoS (Stixrude and Lithgow-Bertelloni, 2005)"
)
eos_names["slb3"] = (
    "Stixrude and Lithgow-Bertelloni 3rd-order EoS (Stixrude and Lithgow-Bertelloni, 2005)"
)
eos_names["slb3-conductive"] = (
    "Stixrude and Lithgow-Bertelloni 3rd-order Conductive EoS (Stixrude and Lithgow-Bertelloni, 2024)"
)
eos_names["dks_l"] = "De Koker and Stixrude Liquid EoS (de Koker and Stixrude, 2013)"
eos_names["dks_s"] = "De Koker and Stixrude Solid EoS (de Koker and Stixrude, 2013)"
eos_names["modular_mgd"] = "Modular Mie-Gruneisen-Debye EoS"
eos_names["modular_mgd_with_anharmonicity"] = (
    "Modular Mie-Gruneisen-Debye EoS with Anharmonicity"
)
eos_names["mt"] = "Modified Tait EoS (Huang and Chow, 1974)"
eos_names["macaw"] = "MACAW EoS (Lozano and Aslam, 2022)"
eos_names["spock"] = "SPOCK EoS (Myhill, 2025)"
eos_names["hp98"] = "Holland and Powell 1998 EoS (Holland and Powell, 1998)"
eos_names["hp_tmt"] = (
    "Holland and Powell Modified Tait Thermal EoS (Holland and Powell, 2011)"
)
eos_names["hp_tmtL"] = (
    "Holland and Powell Modified Tait Thermal Liquid EoS (Holland and Powell, 2011)"
)
eos_names["cork"] = "CORK EoS (Holland and Powell, 1991)"
eos_names["brosh_calphad"] = "Brosh CALPHAD EoS (Brosh, 2007)"
eos_names["aa"] = "Liquid Metal EoS (Anderson and Ahrens, 1994)"
eos_methods = OrderedDict()
eos_methods["murnaghan"] = murnaghan.Murnaghan
eos_methods["bm3shear2"] = bm.BM3Shear2
eos_methods["bm3"] = bm.BM3
eos_methods["bm4"] = bm.BM4
eos_methods["vinet"] = vinet.Vinet
eos_methods["morse"] = morse_potential.Morse
eos_methods["rkprime"] = reciprocal_kprime.RKprime
eos_methods["mgd2"] = mgd.MGD2
eos_methods["mgd3"] = mgd.MGD3
eos_methods["slb2"] = slb.SLB2
eos_methods["slb3"] = slb.SLB3
eos_methods["slb3-conductive"] = slb.SLB3Conductive
eos_methods["dks_l"] = dks_liquid.DKS_L
eos_methods["dks_s"] = dks_solid.DKS_S
eos_methods["modular_mgd"] = mmgd.ModularMGD
eos_methods["modular_mgd_with_anharmonicity"] = mmgd.ModularMGDWithAnharmonicity
eos_methods["mt"] = mt.MT
eos_methods["macaw"] = macaw.MACAW
eos_methods["spock"] = spock.SPOCK
eos_methods["hp98"] = hp.HP98
eos_methods["hp_tmt"] = hp.HP_TMT
eos_methods["hp_tmtL"] = hp.HP_TMTL
eos_methods["cork"] = cork.CORK
eos_methods["brosh_calphad"] = brosh_calphad.BroshCalphad
eos_methods["aa"] = aa.AA
eos_methods["combined"] = CombinedMineralMethod


def create(method):
    """
    Creates an instance of an EquationOfState from a string,
    a class EquationOfState, or an instance of EquationOfState.
    """
    if isinstance(method, str):
        if method in eos_methods:
            return eos_methods[method]()
        else:
            raise Exception("unsupported material method " + method)
    elif isinstance(method, EquationOfState):
        return method
    elif inspect.isclass(method) and issubclass(method, EquationOfState):
        return method()
    else:
        raise Exception("unsupported material method " + method.__class__.__name__)

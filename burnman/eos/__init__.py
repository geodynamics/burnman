# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Equation of State
-----------------


"""


from .equation_of_state import EquationOfState, IsothermalEquationOfState
from .murnaghan import Murnaghan
from .birch_murnaghan import BM3Shear2, BM3, BM4
from .mie_grueneisen_debye import MGD2, MGD3
from .slb import SLB2, SLB3, SLB3Conductive
from .modular_mie_grueneisen_debye import ModularMGD
from .debye_temperature_models import (
    DebyeTemperatureModelBase,
    SLB,
    PowerLawGammaSimple,
    PowerLawGamma,
)
from .modified_tait import MT
from .hp import HP_TMT
from .hp import HP_TMTL
from .hp import HP98
from .cork import CORK
from .vinet import Vinet
from .morse_potential import Morse
from .reciprocal_kprime import RKprime
from .macaw import MACAW
from .spock import SPOCK
from .dks_liquid import DKS_L
from .dks_solid import DKS_S
from .aa import AA
from .brosh_calphad import BroshCalphad

from .property_modifiers import calculate_property_modifications

from .helper import create

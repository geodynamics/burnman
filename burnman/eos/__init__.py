# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Equation of State
-----------------


"""
from __future__ import absolute_import

from .equation_of_state import EquationOfState
from .murnaghan import Murnaghan
from .birch_murnaghan import BM2, BM3
from .birch_murnaghan_4th import BM4
from .mie_grueneisen_debye import MGD2, MGD3
from .slb import SLB2, SLB3
from .modified_tait import MT
from .hp import HP_TMT
from .hp import HP_TMTL
from .hp import HP98
from .cork import CORK
from .vinet import Vinet
from .morse_potential import Morse
from .reciprocal_kprime import RKprime
from .dks_liquid import DKS_L
from .dks_solid import DKS_S
from .aa import AA
from .brosh_calphad import BroshCalphad

from .property_modifiers import calculate_property_modifications

from .helper import create

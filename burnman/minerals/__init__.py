# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
This package contains thermodynamic models from a variety of published sources:
"""


# Stixrude and Lithgow-Bertelloni
from . import SLB_2024
from . import SLB_2022
from . import SLB_2011
from . import SLB_2011_ZSB_2013
from . import SLB_2005

# ab initio
from . import RS_2014_liquids
from . import DKS_2013_liquids
from . import DKS_2013_solids

# Murakami and coworkers
from . import Murakami_etal_2012
from . import Murakami_2013

# Matas and coworkers
from . import Matas_etal_2007

# Holland, Powell and coworkers
from . import HP_2011_ds62
from . import HP_2011_fluids
from . import HHPH_2013
from . import JH_2015
from . import HGP_2018_ds633
from . import ig50NCKFMASHTOCr
from . import ig50NCKFMASTOCr
from . import mb50NCKFMASHTO
from . import mp50KFMASH
from . import mp50MnNCKFMASHTO
from . import mp50NCKFMASHTO

# Kurnosov et al. 2017
from . import KMFBZ_2017

# Irving et al. 2018
from . import ICL_2018

# Other
from . import Sundman_1991
from . import SE_2015
from . import other

# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_user_input_material
---------------------------

Shows user how to input a mineral of his/her choice without usint the library and which physical values
need to be input for BurnMan to calculate :math:`V_P, V_\Phi, V_S` and density at depth.

*Specifically uses:*


* :class:`burnman.mineral.Mineral`

*Demonstrates:*

* how to create your own minerals

"""
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman

# A note about units: all the material parameters are expected to be in plain SI units.
# This means that the elastic moduli should be in Pascals and NOT Gigapascals,
# and the Debye temperature should be in K not C.  Additionally, the reference volume
# should be in m^3/(mol molecule) and not in unit cell volume and 'n' should be
# the number of atoms per molecule.  Frequently in the literature the reference volume
# is given in Angstrom^3 per unit cell.  To convert this to m^3/(mol of molecule)
# you should multiply by 10^(-30) * N_a / Z, where N_a is Avogadro's number and Z is the number of
# atoms per unit cell.  You can look up Z in many places, including
# www.mindat.org

if __name__ == "__main__":

    # input variables ###
    #

    # INPUT for method
    """ choose 'slb2' (finite-strain 2nd order shear modulus,
        stixrude and lithgow-bertelloni, 2005)
    or 'slb3 (finite-strain 3rd order shear modulus,
        stixrude and lithgow-bertelloni, 2005)
    or 'mgd3' (mie-gruneisen-debeye 3rd order shear modulus,
        matas et al. 2007)
    or 'mgd2' (mie-gruneisen-debeye 2nd order shear modulus,
        matas et al. 2007)
    or 'bm2' (birch-murnaghan 2nd order, if you choose to ignore temperature
       (your choice in geotherm will not matter in this case))
    or 'bm3' (birch-murnaghan 3rd order, if you choose to ignore temperature
        (your choice in geotherm will not matter in this case)
    or 'vinet' (vinet equation of state, if you choose to ignore temperature
        (your choice in geotherm will not matter in this case)))"""
    method = 'slb3'

    # in form name_of_mineral (burnman.mineral <- creates list with parameters)
    class own_material (burnman.Mineral):

        def __init__(self):
            self.params = {
                'name': 'myownmineral',
                'equation_of_state': method,
                'V_0': 10.844e-6,  # Molar volume [m^3/(mole molecules)]
                # at room pressure/temperature
                'K_0': 135.19e9,  # Reference bulk modulus [Pa]
                # at room pressure/temperature
                'Kprime_0': 6.04,  # pressure derivative of bulk modulus
                'G_0': 175.0e9,  # reference shear modulus
                # at room pressure/temperature
                'Gprime_0': 1.7,  # pressure derivative of shear modulus
                'molar_mass': .055845,  # molar mass in units of [kg/mol]
                'n': 1,  # number of atoms per formula unit
                'Debye_0': 998.85,  # Debye temperature for material.
                # See Stixrude & Lithgow-Bertelloni, 2005 for values
                'grueneisen_0': 1.368,  # Gruneisen parameter for material.
                # See Stixrude & Lithgow-Bertelloni, 2005 for values
                'q_0': 0.917,  # isotropic strain derivative of gruneisen
                # parameter. Values in Stixrude & Lithgow-Bertelloni, 2005
                'eta_s_0': 3.0  # full strain derivative of gruneisen parameter
                # parameter. Values in Stixrude & Lithgow-Bertelloni, 2005
            }
            burnman.Mineral.__init__(self)

    rock = own_material()

    # seismic model for comparison: (see burnman/seismic.py)
    seismic_model = burnman.seismic.PREM()  # pick from .prem() .slow() .fast()
    number_of_points = 20  # set on how many depth slices the computations should be done
    depths = np.linspace(700e3, 2800e3, number_of_points)
    # depths = seismic_model.internal_depth_list(mindepth=700.e3,
    # maxdepth=2800.e3)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ['pressure', 'density', 'v_p', 'v_s', 'v_phi'], depths)

    temperature = burnman.geotherm.brown_shankland(seis_p)
    # The next line is not required here, because the method is set
    # automatically by defining 'equation_of_state' in mineral.params. This
    # shows an alternative way to set the method later, or reset the method to
    # a different one.
    rock.set_method(method)

    print("Calculations are done for:")
    rock.debug_print()

    mat_rho, mat_vs, mat_vphi = \
        rock.evaluate(['density', 'v_s', 'v_phi'], seis_p, temperature)

    [vs_err, vphi_err, rho_err] = \
        burnman.compare_chifactor(
            [mat_vs, mat_vphi, mat_rho], [seis_vs, seis_vphi, seis_rho])

    print(vs_err, vphi_err, rho_err)

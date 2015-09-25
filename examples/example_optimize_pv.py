# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.


"""
    
example_optimize_pv
-------------------

Vary the amount perovskite vs. ferropericlase and compute the error in the
seismic data against PREM. For more extensive comments on this setup, see tutorial/step_2.py

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :class:`burnman.seismic.PREM`
* :func:`burnman.geotherm.brown_shankland`
* :func:`burnman.main.velocities_from_rock`
* :func:`burnman.main.compare_l2`

*Demonstrates:*

* compare errors between models
* loops over models

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":


### input variables ###
    #######################


    seismic_model = burnman.seismic.PREM() # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    number_of_points = 20 #set on how many depth slices the computations should be done
    depths = np.linspace(700e3,2800e3, number_of_points)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(['pressure','density','v_p','v_s','v_phi'],depths)

    temperature = burnman.geotherm.brown_shankland(seis_p)

    def material_error(amount_perovskite):
        #Define composition using the values from Murakami et al. 2012 (Note: fe_perovskite and fe_periclase do not represent pure iron
        #endmembers here, but contain 6% and 20% Fe respectively. 
        rock = burnman.Composite([amount_perovskite, 1.0-amount_perovskite],
                                 [minerals.Murakami_etal_2012.fe_perovskite(),
                                  minerals.Murakami_etal_2012.fe_periclase()])


        mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = \
            burnman.velocities_from_rock(rock, seis_p, temperature, burnman.averaging_schemes.VoigtReussHill())

        print "Calculations are done for:"
        rock.debug_print()

        [vs_err, vphi_err, rho_err] = \
            burnman.compare_l2(depths, [mat_vs,mat_vphi,mat_rho], [seis_vs,seis_vphi,seis_rho])

        return vs_err, vphi_err

    xx=np.linspace(0.0, 1.0, 40)
    errs=np.array([material_error(x) for x in xx])
    yy_vs=errs[:,0]
    yy_vphi=errs[:,1]
    plt.plot (xx,yy_vs,"r-x",label=("vs error"))
    plt.plot (xx,yy_vphi,"b-x",label=("vphi error"))
    plt.yscale('log')
    plt.xlabel('% Perovskite')
    plt.ylabel('Error')
    plt.legend()
    plt.show()

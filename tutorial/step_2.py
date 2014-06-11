# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
"""

import numpy as np
import matplotlib.pyplot as plt

import burnman
from burnman import minerals


n_depths = 20
min_depth = 850.e3
max_depth = 2800.e3
depths = np.linspace(min_depth, max_depth, n_depths)

seismic_model = burnman.seismic.PREM()
pressure, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)
temperature = burnman.geotherm.brown_shankland(pressure)

def material_error(amount_perovskite):
    rock =  ?????          '''<-------------------- Create a rock'''
    rock.set_method('slb3')
    density, vp, vs, vphi, K, G =   ???? '''<---------------- Calculate the elastic properties of hte rock'''

    print "Calculations are done for:"
    rock.debug_print()

    [vs_err, vphi_err, rho_err]=burnman.compare_l2(depths,[vs,vphi,density],[seis_vs,seis_vphi,seis_rho])

    return vs_err, vphi_err, rho_err

xx=np.linspace(0.0, 1.0, 50)
errs=np.array([material_error(x) for x in xx])
yy_vs=errs[:,0]
yy_vphi=errs[:,1]
yy_rho=errs[:,2]

plt.plot (xx,yy_vs,"r-x",label=("vs error"))
plt.plot (xx,yy_vphi,"b-x",label=("vphi error"))
plt.plot (xx,yy_rho,"g-x",label=("rho error"))
plt.yscale('log')
plt.xlabel('% Perovskite')
plt.ylabel('Error')
plt.legend()
plt.show()

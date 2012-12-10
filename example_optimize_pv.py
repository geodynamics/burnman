import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
	sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

### input variables ###
#######################

#INPUT for method
method = 'slb' # choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)

seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see code/seismic.py)
number_of_points = 20 #set on how many depth slices the computations should be done
depths = np.linspace(700,2800, number_of_points)
seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

geotherm = burnman.geotherm.brown_shankland
temperature = [geotherm(p) for p in seis_p]

def material_error(amount_perovskite):
	phases = [minerals.Murakami_perovskite(), minerals.Murakami_fp_LS()]
	molar_abundances = [amount_perovskite, 1.0-amount_perovskite]

	for ph in phases:
		ph.set_method(method)

	print "Calculations are done for:"
	for i in range(len(phases)):
		print molar_abundances[i], " of phase", phases[i].to_string()

	mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu = burnman.calculate_velocities(seis_p, temperature, phases, molar_abundances)	

	#[rho_err,vphi_err,vs_err]=burnman.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho)
	[rho_err,vphi_err,vs_err]=burnman.compare_two(depths,mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho)

	return vs_err, vphi_err

xx=np.linspace(0.0, 1.0, 40)
errs=np.array([material_error(x) for x in xx])
yy_vs=errs[:,0]
yy_vphi=errs[:,1]
plt.plot (xx,yy_vs,label=("vs error"))
plt.plot (xx,yy_vphi,label=("vphi error"))
plt.yscale('log')

#plt.ylim(0,100)
plt.xlabel('% Perovskite')
plt.ylabel('Error')
plt.legend()
plt.show()

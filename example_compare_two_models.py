import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
	sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

###Input Model 1

##input method
method_1 = 'mgd' # choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)


#Input composition of model 1. See example_composition for potential choices. We'll just choose something simple here
	
phases_1 = (minerals.Murakami_perovskite(), minerals.Murakami_fp_LS())
amount_perovskite_1 = 0.95
molar_abundances_1 = ( amount_perovskite_1, 1.0-amount_perovskite_1)

#input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array

seis_p_1 = np.arange(25e9, 125e9, 5e9)

#input your geotherm. Either choose one (See example_geotherms.py) or create one.We'll use Brown and Shankland.

geotherm = burnman.geotherm.brown_shankland
temperature_1 = [geotherm(p) for p in seis_p_1]

##Now onto the second model parameters

##input second method
method_2 = 'slb' # choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)


#Input composition of model 2. See example_composition for potential choices. We'll just choose something simple here
	
phases_2 = (minerals.Murakami_perovskite(), minerals.Murakami_fp_LS())
amount_perovskite_2 = 0.95
molar_abundances_2 = ( amount_perovskite_2, 1.0-amount_perovskite_2)

#input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array

seis_p_2 = np.arange(25e9, 125e9, 5e9)

#input your geotherm. Either choose one (See example_geotherms.py) or create one.We'll use Brown and Shankland.

geotherm = burnman.geotherm.brown_shankland
temperature_2 = [geotherm(p) for p in seis_p_2]


#Now we'll calculate the models. 

for ph in phases_1:
	ph.set_method(method_1)

print "Calculations are done for:"
for i in range(len(phases_1)):
	print molar_abundances_1[i], " of phase", phases_1[i].to_string()

mat_rho_1, mat_vs_1, mat_vp_1, mat_vphi_1, mat_K_1, mat_mu_1 = burnman.calculate_velocities(seis_p_1, temperature_1, phases_1, molar_abundances_1)	

for ph in phases_2:
	ph.set_method(method_2)

print "Calculations are done for:"
for i in range(len(phases_2)):
	print molar_abundances_2[i], " of phase", phases_2[i].to_string()

mat_rho_2, mat_vs_2, mat_vp_2, mat_vphi_2, mat_K_2, mat_mu_2 = burnman.calculate_velocities(seis_p_2, temperature_2, phases_2, molar_abundances_2)	


##Now let's plot the comparison. You can conversely just output to a data file (see example_woutput.py)

plt.subplot(2,2,1)
p1,=plt.plot(seis_p_1/1.e9,mat_vs_1,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=4)
p2,=plt.plot(seis_p_2/1.e9,mat_vs_2,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
plt.title("Vs (km/s)")


# plot Vphi
plt.subplot(2,2,2)
p1,=plt.plot(seis_p_1/1.e9,mat_vphi_1,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=4)
p2,=plt.plot(seis_p_2/1.e9,mat_vphi_2,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
plt.title("Vphi (km/s)")

# plot density
plt.subplot(2,2,3)
p1,=plt.plot(seis_p_1/1.e9,mat_rho_1,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=4)
p2,=plt.plot(seis_p_2/1.e9,mat_rho_2,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
plt.title("density (kg/m^3)")
plt.legend([p1,p2],["model 1", "model 2"], loc=4)


#plot percent differences
mat_vs_2_interp = np.interp(seis_p_1, seis_p_2,mat_vs_2)
mat_vphi_2_interp = np.interp(seis_p_1, seis_p_2,mat_vphi_2)
mat_rho_2_interp = np.interp(seis_p_1, seis_p_2,mat_rho_2)

per_diff_vs = 100*(mat_vs_1 - mat_vs_2_interp)/mat_vs_1
per_diff_vphi = 100*(mat_vphi_1 - mat_vphi_2_interp)/mat_rho_1
per_diff_rho = 100*(mat_rho_1 - mat_rho_2_interp)/mat_rho_1

plt.subplot(2,2,4)
p1,=plt.plot(seis_p_1/1.e9,per_diff_vs,color='g',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
p2,=plt.plot(seis_p_1/1.e9,per_diff_vphi,color='c',linestyle='-',marker='+',markerfacecolor='k',markersize=4)
p3,=plt.plot(seis_p_1/1.e9,per_diff_rho,color='m',linestyle='-',marker='x',markerfacecolor='k',markersize=4)

plt.title("percent difference")
plt.legend([p1,p2,p3],["vs", "vphi","density"], loc=4)

plt.show()

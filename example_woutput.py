import os, sys, numpy as np
import matplotlib.pyplot as plt
import csv

lib_path = os.path.abspath('code/')
sys.path.append(lib_path)

from code import minerals 
from code import main as main

#See bottom for output tool


### input variables ###
#######################

#INPUT for method
method = 'slb' # choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)

#specify material
phases = (minerals.mg_fe_perovskite(0.7), minerals.ferropericlase(0.5))
amount_perovskite = 0.95
molar_abundances = ( amount_perovskite, 1.0-amount_perovskite)

#define some pressure range
pressures = np.arange(25e9,130e9,5e9)

geotherm = main.get_geotherm("brown_shankland")
temperature = [geotherm(p) for p in pressures]

for ph in phases:
	ph.set_method(method)

print "Calculations are done for:"
for i in range(len(phases)):
	print molar_abundances[i], " of phase", phases[i].to_string()

mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu = main.calculate_velocities(pressures, temperature, phases, molar_abundances)	
	
#write to file:
output_filename = "example_woutput.txt" 
file = open(output_filename, 'wb')
file.write("#Pressure\tTemperature\tmat_rho\tmat_vs\tmat_vp\tmat_vphi\tmat_K\tmat_mu\n")

data = zip(pressures,temperature,mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu)
np.savetxt(file,data,fmt='%.10e',delimiter='\t')

	
print "\nYour data has been saved as: ",output_filename


#standard numpy, scipy imports
import numpy as np
import matplotlib.pyplot as plt


#eos imports
import seismo_in as seis
from  minerals import *
import comparison as comp
import geotherm as gt
import composition as part

# read in userinput
import sys
import imp
file=sys.argv[1]+"/userinput.py"
print file
userinput=(imp.load_source("userinput",file))


def plot_eos(phases, molar_abundances, geotherm):
        # read in seismic model
	s=seis.seismo_in()
        print "Looking up seismic model:",s.params['name']
        seis_p = np.array(s.pressure())*1.e9     
        seis_r = np.array(s.radius())/1.e3
        seis_vp = np.array(s.v_p())/1.e3#[seis.seis_V(y/1.e9)[0]/1.e3 for y in seis_p]
        seis_vphi = np.array(s.v_phi())/1.e3#
        seis_vs = np.array(s.v_s())/1.e3#[seis.seis_V(y/1.e9)[1]/1.e3 for y in seis_p]
        seis_rho = np.array(s.density())#[seis.seis_density(y/1.e9) for y in seis_p]
 
        
	temperature = [geotherm(p) for p in seis_p]
	mat_vs = np.empty_like(seis_p)
	mat_vp = np.empty_like(seis_p)
        mat_vphi = np.empty_like(seis_p)
	mat_rho = np.empty_like(seis_p)
	mat_K = np.empty_like(seis_p)
	mat_mu = np.empty_like(seis_p)

	print "Calculating elastic properties for phase assemblage \n"
	print "seismic p (GPa)	T (K)	density(kg/m^3)	K(Gpa) G(GPa)	Vs (km/s)	Vp(km/s)	Vphi (km/s)"
	for i in range(len(seis_p)):
		rho,vp,vs,vphi,K,mu = \
		vrh.voigt_reuss_hill(seis_p[i], temperature[i], phases, molar_abundances)
		print seis_p[i]/1.e9,"	", temperature[i],"	",rho,"	", K,"	", mu,"	", vs/1.e3,"	", vp/1.e3,"	", vphi/1.e3
		#print seis_p[i]/1.e9,"	",rho,"	", mu
		mat_rho[i] = rho/1.e3
		mat_vs[i] = vs/1.e3
		mat_vp[i] = vp/1.e3
                mat_vphi[i] = vphi/1.e3
                mat_K[i] = K
                mat_mu[i] = mu

	[rho_err,vp_err,vs_err]=comp.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho)
            
	plt.subplot(2,2,1)
	p1,=plt.plot(seis_p/1.e9,mat_vs,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4)
	p2,=plt.plot(seis_p/1.e9,seis_vs,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
	plt.title("Vs (km/s)")
	plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
	plt.ylim(6.,7.5)
        plt.text(40,7.3,"misfit= %3.3f" % vs_err)

	# plot Vphi
	plt.subplot(2,2,2)
	p1,=plt.plot(seis_p/1.e9,mat_vphi,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4)
	p2,=plt.plot(seis_p/1.e9,seis_vphi,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
	plt.title("Vphi (km/s)")
	plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
	plt.ylim(7,12)
#	plt.legend([p1,p2],["Murakami (0.93 Pv, 0.07 fp)", "seismic model (PREM)"], loc=4)
	plt.text(40,14.5,"misfit= %3.3f" % vp_err)

	# plot density
	plt.subplot(2,2,3)
	p1,=plt.plot(seis_p/1.e9,mat_rho,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4)
	p2,=plt.plot(seis_p/1.e9,seis_rho/1.e3,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
	plt.title("density (kg/m^3)")
	plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
        plt.text(40,5.3,"misfit= %3.3f" % rho_err)

        # plot geotherm
	plt.subplot(2,2,4)
	p1,=plt.plot(seis_p/1e9,temperature,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=4)
	plt.title("Geotherm (K)")
	plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
	


        plt.savefig(sys.argv[1]+"/output_plots.png")
	plt.show()

##### main


#import sys
#import imp
#userinput=imp.load_source("userinput",sys.argv[1])

# read in geotherm
if (userinput.geotherm=='geotherm_brown_shankland'):
	geothermal = lambda p: gt.geotherm_brown_shankland(p)
	print "geotherm is from Brown & Shankland (1981)"
else:
	print "define geothermal gradient"
# read in composition
if (userinput.composition_input=='weight_percents'):
	weight_percents = userinput.weight_percents 
	phase_fractions,relative_molar_percent= part.conv_inputs(weight_percents)
	iron_content = lambda p,t: part.calculate_partition_coefficient(p,t,relative_molar_percent)
	pv_d = mg_fe_perovskite_pt_dependent(iron_content)
	fp_d = ferropericlase_pt_dependent(iron_content)
        phases = (pv_d, fp_d)
	iron_content=("(P,T) dependent","(P,T) dependent")	
	molar_abundances = (phase_fractions['pv'],phase_fractions['fp'])
elif (userinput.composition_input=='2phase_fractions'):
	if(userinput.calculate_partitioning=='auto'):
		iron_content = lambda p,t: part.calculate_partition_coefficient(p,t,relative_molar_percent)
		pv_d = mg_fe_perovskite_pt_dependent(iron_content)
		fp_d = ferropericlase_pt_dependent(iron_content)
        	pv = eval(userinput.pv)(pv_d)
        	fp = eval(userinput.fp)(fp_d)
        elif (userinput.calculate_partitioning=='on'):
		pv_d=userinput.partitioning['pv']
		fp_d=userinput.partitioning['fp']
		iron_content=(pv_d,fp_d)
        	pv = eval(userinput.pv)(pv_d)
        	fp = eval(userinput.fp)(fp_d)
	else:
		pv = eval(userinput.pv)()
  		fp = eval(userinput.fp)()
		iron_content=('not used', 'not used')
        phases = (pv, fp)
	molar_abundances = (userinput.phase_fractions['pv'],userinput.phase_fractions['fp'])
elif (userinput.composition_input=='nphase_fractions'):
        list=userinput.phase_names
        molar_abundances=np.empty(len(list))
        iron_content=np.empty(len(list))
        phases=[]
	for i in range(len(list)):
		try:
			molar_abundances[i]=float(userinput.phase_fractions[list[i]]) 
		except:
			molar_abundances[i]=0.
		try:	
			iron_content[i]=userinput.partitioning[list[i]]
		except:
			iron_content[i]=None
		print iron_content[i]
		if np.isnan(iron_content[i]):
			try: 
			        phase=eval(str(userinput.phases[i]))()
			        print phase
			except:
				raise Exception("Iron partitioning needs to be defined for", userinput.phases[i])
		else:
			try:
                		phase=eval(str(userinput.phases[i]))(iron_content[i])
                	except:
				raise Exception(userinput.phases[i], "does not handle an iron_partioning #")		
								
		phases.append(phase)
else:
	raise("Choose method to determine composition: 'weight_percents' or 'phase_fractions'")

print "Calculations are done for:"
for i in range(len(phases)):
	print molar_abundances[i], " of phase", phases[i], " with partitioning fraction", iron_content[i]

plot_eos(phases, molar_abundances, geothermal)


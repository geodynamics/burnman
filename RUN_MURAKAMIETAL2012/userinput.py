### input variables ###
#######################

#INPUT for method
method = 'slb'					# choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005, ONLY ONE WORKING) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)


#INPUT for geotherm
geotherm = 'geotherm_brown_shankland'



#INPUT for composition
# for material constants see minerals.py
pv = 'Murakami_perovskite'            	# choose 'mg_fe_perovskite' or 'Murakami_perovskite'
fp = 'Murakami_fp_LS'                   	# choose 'ferropericlase', or 'Murakami_fp_LS'

# for material constants see minerals.py
composition_input = '2phase_fractions'		# choose 'weight_percents', '2phase_fractions', or 'nphase_fractions'

# if '2phase_fractions' or 'nphase_fractions'
phase_names=('pv','fp')
phases = (pv,fp)
phase_fractions = {'pv':0.93, 'fp':0.07} 	# should add up to 1.0

# if '2phase_fractions'
calculate_partitioning = 'off'          	# sets if partioning coeefficients are used or calculated, 'on'/'off'/'auto', only for 2 phase fractions
partitioning = {'pv':0.94, 'fp':0.79}   	# ignored if calculate_partioning = 'on'artitioning = {'pv':0.94, 'fp':0.79}   # ignored if calculate_partioning = 'on'

#INPUT for seismic models
name = 'prem'  					# choose from 'prem'(at 1 second, Dziewonski & Anderson, 1981), 'ref_fast' (Lekic et al. 2012), 'ref_slow' (Lekic et al. 2012)
attenuation_correction= 'off'   		# correction for attuation (not required for PREM at 1s, see Matas et al. 2007, page 4) choose from 'on', 'off'
depth_min = 0.0   				# minimum depth considered in km, choose 0.0 to use the limits given by the seismic model
depth_max = 0.0   				# maximum depth considered in km, choose 0.0 to use the limits given by the seismic model
depth_step = 100.0 				# steps at which the seismic velocity and calculated velocity are compared (in between is interpolated)


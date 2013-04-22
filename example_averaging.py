# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example shows the effect of different averaging schemes.

requires:
- geotherms
- compute seismic velocities

teaches:
- averaging

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

if __name__ == "__main__":    
    
    #INPUT for method
    """ choose 'slb2' (finite-strain 2nd order sheer modulus, stixrude and lithgow-bertelloni, 2005)
    or 'slb3 (finite-strain 3rd order shear modulus, stixrude and lithgow-bertelloni, 2005)
    or 'mgd3' (mie-gruneisen-debeye 3rd order shear modulus, matas et al. 2007)
    or 'mgd2' (mie-gruneisen-debeye 2nd order shearl modulus, matas et al. 2007)
    or 'bm' (birch-murnaghan, if you choose to ignore temperature (your choice in geotherm will not matter in this case))"""
    
    method = 'slb3'

#     class olivine (burnman.material):
#         """
#         Stixrude & Lithgow-Bertelloni 2005 and references therein 
#         """
#         def __init__(self):
#             self.params = {
#             'equation_of_state':'slb3',
#             'ref_V': 46.29e-6,
#             'ref_K': 135.0e9,
#             'K_prime': 4.2,
#             'ref_mu': 51.0e9,
#             'mu_prime': 1.4,
#             'molar_mass': .1406,
#             'n': 6,
#             'ref_Debye': 619.,
#             'ref_grueneisen': 1.06,
#             'q0': 3.6,
#             'eta_0s': 1.1 }     
    
    amount_perovskite = 0.6
    #rock = burnman.composite( ( (minerals.mg_perovskite(), amount_perovskite), 
    #                                (olivine(), 1.0-amount_perovskite) ) )
    rock = burnman.composite( ( (minerals.mg_perovskite(), amount_perovskite), 
                                    (minerals.periclase(), 1.0-amount_perovskite) ) )
    #rock = burnman.composite ( ( (minerals.Murakami_fe_perovskite(), amount_perovskite),
    #                                 (minerals.Murakami_fe_periclase_LS(), 1.0-amount_perovskite)) )
   
    #rock = burnman.composite ( ( (minerals.mg_perovskite(), amount_perovskite),
    #                                 (minerals.Matas_periclase(), 1.0-amount_perovskite)) )
           
    #seismic model for comparison:
    seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    number_of_points = 20 #set on how many depth slices the computations should be done
    # we will do our computation and comparison at the following depth values:
    depths = np.linspace(700e3, 2800e3, number_of_points)
    #alternatively, we could use the values where prem is defined:
    #depths = seismic_model.internal_depth_list()
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)
    
            
    geotherm = burnman.geotherm.brown_shankland
    temperature = [geotherm(p) for p in seis_p]
    
    rock.set_method(method)
    
    print "Calculations are done for:"
    for ph in rock.phases:
        print ph.fraction, " of phase", ph.mineral.to_string()
    
    moduli_list = burnman.calculate_moduli(rock, seis_p, temperature)

    moduli = burnman.average_moduli(moduli_list, burnman.averaging_schemes.voigt_reuss_hill())
    mat_vs1, mat_vp1, mat_vphi1 = burnman.compute_velocities(moduli)
    mat_K1, mat_mu1, mat_rho1 = moduli.K, moduli.mu, moduli.rho

    moduli = burnman.average_moduli(moduli_list, burnman.averaging_schemes.linear())
    mat_vs2, mat_vp2, mat_vphi2 = burnman.compute_velocities(moduli)
    mat_K2, mat_mu2, mat_rho2 = moduli.K, moduli.mu, moduli.rho

    moduli = burnman.average_moduli(moduli_list, burnman.averaging_schemes.voigt())
    mat_vs3, mat_vp3, mat_vphi3 = burnman.compute_velocities(moduli)
    mat_K3, mat_mu3, mat_rho3 = moduli.K, moduli.mu, moduli.rho

    moduli = burnman.average_moduli(moduli_list, burnman.averaging_schemes.reuss())
    mat_vs4, mat_vp4, mat_vphi4 = burnman.compute_velocities(moduli)
    mat_K4, mat_mu4, mat_rho4 = moduli.K, moduli.mu, moduli.rho


    moduli = moduli_list[0]
    mat_vsa, mat_vpa, mat_vphia = burnman.compute_velocities(moduli)
    mat_Ka, mat_mua, mat_rhoa = moduli.K, moduli.mu, moduli.rho

    moduli = moduli_list[1]
    mat_vsb, mat_vpb, mat_vphib = burnman.compute_velocities(moduli)
    mat_Kb, mat_mub, mat_rhob = moduli.K, moduli.mu, moduli.rho
    
    #[rho_err,vphi_err,vs_err]=burnman.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho)
        
    
    # PLOTTING
    
    # plot vs
    plt.plot(seis_p/1.e9,mat_vs3/1.e3,color='c',linestyle='-',marker='x',markersize=4,label='voigt')
    plt.plot(seis_p/1.e9,mat_vs4/1.e3,color='k',linestyle='-',marker='x',markersize=4,label='reuss')
    plt.plot(seis_p/1.e9,mat_vs1/1.e3,color='b',linestyle='-',marker='x',markersize=4,label='VRH')
    plt.plot(seis_p/1.e9,mat_vs2/1.e3,color='r',linestyle='-',marker='x',markersize=4,label='linear')
    plt.plot(seis_p/1.e9,mat_vsa/1.e3,color='y',linestyle='--',marker='x',markersize=4,label='mg perovskite')
    plt.plot(seis_p/1.e9,mat_vsb/1.e3,color='g',linestyle='--',marker='x',markersize=4,label='periclase')
    plt.title("Vs (km/s), mix: 60% mg perovskite")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    #plt.ylim(5.2,7.0)
    plt.legend(loc='lower right')

    #plt.savefig("output_figures/example_averaging.png")
    plt.show()

    mat_vsa_norm=(mat_vsa-mat_vsb)/(mat_vsa-mat_vsb)
    mat_vsb_norm=(mat_vsb-mat_vsb)/(mat_vsa-mat_vsb)
    mat_vs1_norm=(mat_vs1-mat_vsb)/(mat_vsa-mat_vsb)
    mat_vs2_norm=(mat_vs2-mat_vsb)/(mat_vsa-mat_vsb)
    mat_vs3_norm=(mat_vs3-mat_vsb)/(mat_vsa-mat_vsb)
    mat_vs4_norm=(mat_vs4-mat_vsb)/(mat_vsa-mat_vsb)
    plt.plot(seis_p/1.e9,mat_vs3_norm,color='c',linestyle='-',marker='x',markersize=4,label='voigt')
    plt.plot(seis_p/1.e9,mat_vs4_norm,color='k',linestyle='-',marker='x',markersize=4,label='reuss')
    plt.plot(seis_p/1.e9,mat_vs1_norm,color='b',linestyle='-',marker='x',markersize=4,label='VRH')
    plt.plot(seis_p/1.e9,mat_vs2_norm,color='r',linestyle='-',marker='x',markersize=4,label='linear')
    plt.plot(seis_p/1.e9,mat_vsa_norm,color='y',linestyle='--',marker='x',markersize=4,label='mg perovskite')
    plt.plot(seis_p/1.e9,mat_vsb_norm,color='g',linestyle='--',marker='x',markersize=4,label='periclase')
    plt.title("Vs (km/s), mix: 60% mg perovskite")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    #plt.ylim(5.2,7.0)
    plt.legend(loc='lower right')

    #plt.savefig("output_figures/example_averaging_normalized.png")
    plt.show()

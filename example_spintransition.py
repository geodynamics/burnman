# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example shows the different minerals that are implemented with a spin
transition.  Minerals with spin transition are implemented by defining two
separate minerals (one for the low and one for the high spin state).  Then a
third dynamic mineral is created that switches between the two previously
defined minerals by comparing the current pressure to the transition pressure.

requires:
- geotherms
- defined minerals, from library or user created (see: example_user_input_material.py)
- seismic models
- compute seismic velocities

teaches:
- implementation of spin transition in (Mg,Fe)O at user defined pressure
"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

if __name__ == "__main__":    
    
    #seismic model for comparison:
    # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    seismic_model = burnman.seismic.prem() 
    number_of_points = 40 #set on how many depth slices the computations should be done
    # we will do our computation and comparison at the following depth values:
    depths = np.linspace(700e3, 2800e3, number_of_points)
    #alternatively, we could use the values where prem is defined:
    #depths = seismic_model.internal_depth_list()
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

    # here we use the Brown & Shankland geotherm    
    geotherm = burnman.geotherm.brown_shankland
    temperature = [geotherm(p) for p in seis_p]    
    
    
    # We create one mineral that contains spin transitions
    rock = burnman.composite ( [(minerals.Muretal2012.fe_periclase(), 1.0)] )

    # The mineral Murakami_fe_periclase is derived from minerals.helper_spin_transition
    # which contains the logic to switch between two other minerals based on the
    # current pressure. The mineral is implemented similar to the following lines:
    #
    #   class Murakami_fe_periclase(helper_spin_transition):
    #     def __init__(self):
    #       helper_spin_transition.__init__(self, 63.0e9, Murakami_fe_periclase_LS(), Murakami_fe_periclase_HS())
    #
    # Note the reference to the low spin and high spin minerals (_LS and _HS).
    
    # Set method, here set to 'slb2' as the shear wave moduli in
    # Murakami et al. 2012 were fit to second order 
    rock.set_method('slb2')
    
    # Now we calculate the velocities
    print "Calculations are done for:"
    for ph in rock.phases:
        print ph.fraction, " of phase", ph.mineral.to_string()
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_mu = \
       burnman.velocities_from_rock(rock, seis_p, temperature, burnman.averaging_schemes.voigt_reuss_hill())
       
    # plot example 1
    plt.subplot(2,2,1)
    plt.plot(seis_p/1.e9,mat_vs/1.e3,color='b',linestyle='-',marker='o',\
    markerfacecolor='b',markersize=4,label='Vs')
    plt.plot(seis_p/1.e9,mat_vphi/1.e3,color='r',linestyle='-',marker='o',\
    markerfacecolor='r',markersize=4, label='Vp')
    plt.plot(seis_p/1.e9,mat_rho/1.e3,color='k',linestyle='-',marker='o',\
    markerfacecolor='k',markersize=4, label='rho')
    plt.title("ferropericlase (Murakami et al. 2012)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.ylim(5,12)
    plt.legend(loc='upper left')
    
    # example 2: Here we show the effects of using purely High Spin or Low Spin
    
    rock = burnman.composite( [(minerals.Muretal2012.fe_periclase_LS(), 1.0)] )
    rock.set_method('slb2')
    
    mat_rho_LS, mat_vp_LS, mat_vs_LS, mat_vphi_LS, _, _ = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
        burnman.averaging_schemes.voigt_reuss_hill())
    
    rock = burnman.composite( [(minerals.Muretal2012.fe_periclase_HS(), 1.0)] )
    rock.set_method('slb2')
    mat_rho_HS, mat_vp_HS, mat_vs_HS, mat_vphi_HS, _, _ = \
            burnman.velocities_from_rock(rock, seis_p, temperature, \
            burnman.averaging_schemes.voigt_reuss_hill())

    
    rock = burnman.composite( [(minerals.Muretal2012.fe_periclase(), 1.0)] )
    rock.set_method('slb2')
    mat_rho_ON, mat_vp_ON, mat_vs_ON, mat_vphi_ON, _, _ = \
            burnman.velocities_from_rock(rock, seis_p, temperature, \
            burnman.averaging_schemes.voigt_reuss_hill())
    
    plt.subplot(2,2,2)
    plt.plot(seis_p/1.e9,mat_vs_LS/1.e3,color='b',linestyle='-',marker='.',\
    markerfacecolor='b',markersize=4,label='Vs LS')
    plt.plot(seis_p/1.e9,mat_vs_HS/1.e3,color='r',linestyle='-',marker='.',\
    markerfacecolor='b',markersize=4,label='Vs HS')
    plt.plot(seis_p/1.e9,mat_vs_ON/1.e3,color='g',linestyle='-',marker='o',\
    markerfacecolor='b',markersize=4,label='Vs ON')
    plt.title("Murakami_fp")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    #plt.ylim(300,800)
    plt.legend(loc='lower right')
    
    
    # Example 3: Periclase from Speziale et al. 2006 
    # Here the compositions are implemented as fixed minerals. 
    #For other options see example_composition.py
    rock = burnman.composite( [(minerals.other.Speziale_fe_periclase(), 1.0)] )
    
    rock.set_method('slb3')
    
    print "Calculations are done for:"
    for ph in rock.phases:
        print ph.fraction, " of phase", ph.mineral.to_string()
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_mu = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
        burnman.averaging_schemes.voigt_reuss_hill())
    
    # plot example 3
    plt.subplot(2,2,3)
    plt.plot(seis_p/1.e9,mat_rho/1.e3,color='k',linestyle='-',marker='o',\
    markerfacecolor='k',markersize=4, label='rho')
    plt.title("ferropericlase (Speziale et al. 2007)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.legend(loc='upper left')
    
    plt.subplot(2,2,4)
    plt.plot(seis_p/1.e9,mat_vphi/1.e3,color='b',linestyle='-',marker='o',\
    markerfacecolor='b',markersize=4,label='Vphi')
    plt.title("ferropericlase (Speziale et al. 2007)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.legend(loc='upper left') 
    
    plt.savefig("output_figures/examples_spintransition.png")
    plt.show()

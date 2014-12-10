# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
    
example_compare_all_methods
---------------------------

This example demonstrates how to call each of the individual calculation
methodologies that exist within BurnMan. See below for current options. This
example calculates seismic velocity profiles for the same set of minerals and
a plot of Vs, Vphi and density is produce for the user to compare each of the
different methods.

*Specifically uses:*

* :doc:`eos`


*Demonstrates:*

* Each method for calculating velocity profiles currently included within BurnMan

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":
    #Input composition.

    amount_perovskite = 0.95
    rock = burnman.Composite([amount_perovskite, 1.0-amount_perovskite],
                             [minerals.Murakami_etal_2012.fe_perovskite(),
                              minerals.Murakami_etal_2012.fe_periclase_LS()])

    #(min pressure, max pressure, pressure step)
    seis_p = np.arange(25e9, 125e9, 5e9)

    #Input adiabat potential temperature
    T0 = 1500.0

    #Now we'll calculate the models by forcing the rock to use a method. The preset eqauation of state for the Murakami_etal_2012 minerals is 'slb2'
    
    """ 'slb2' (finite-strain 2nd order shear modulus,
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
        (your choice in geotherm will not matter in this case))"""

    rock.set_method('mgd3')
    temperature = burnman.geotherm.adiabatic(seis_p, T0, rock)

    print "Calculations are done for:"
    rock.debug_print()

    mat_rho_1, mat_vp_1, mat_vs_1, mat_vphi_1, mat_K_1, mat_G_1 = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
                                     burnman.averaging_schemes.VoigtReussHill())

    rock.set_method('slb2')
    temperature = burnman.geotherm.adiabatic(seis_p, T0, rock)

    mat_rho_2, mat_vp_2, mat_vs_2, mat_vphi_2, mat_K_2, mat_G_2 = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
                                     burnman.averaging_schemes.VoigtReussHill())

    rock.set_method('slb3')
    temperature = burnman.geotherm.adiabatic(seis_p, T0, rock)

    mat_rho_3, mat_vp_3, mat_vs_3, mat_vphi_3, mat_K_3, mat_G_3 = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
                                     burnman.averaging_schemes.VoigtReussHill())

    rock.set_method('bm2')
    temperature = burnman.geotherm.adiabatic(seis_p, T0, rock)

    mat_rho_4, mat_vp_4, mat_vs_4, mat_vphi_4, mat_K_4, mat_G_4 = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
                                     burnman.averaging_schemes.VoigtReussHill())

    rock.set_method('bm3')
    temperature = burnman.geotherm.adiabatic(seis_p, T0, rock)
    mat_rho_5, mat_vp_5, mat_vs_5, mat_vphi_5, mat_K_5, mat_G_5 = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
                                     burnman.averaging_schemes.VoigtReussHill())

    ##Now let's plot the comparison. You can conversely just output to a data file
    #(see example_woutput.py)

    plt.subplot(2,2,1)
    plt.plot(seis_p/1.e9,mat_vs_1/1.e3,color='r',linestyle='-',marker='^', \
             markerfacecolor='r',markersize=4)
    plt.plot(seis_p/1.e9,mat_vs_2/1.e3,color='k',linestyle='-',marker='v', \
             markerfacecolor='k',markersize=4)
    plt.plot(seis_p/1.e9,mat_vs_3/1.e3,color='b',linestyle='-',marker='x', \
             markerfacecolor='b',markersize=4)
    plt.plot(seis_p/1.e9,mat_vs_4/1.e3,color='g',linestyle='-',marker='o', \
             markerfacecolor='g',markersize=4)
    plt.plot(seis_p/1.e9,mat_vs_5/1.e3,color='y',linestyle='-',marker='*', \
             markerfacecolor='y',markersize=4)
    plt.title("Vs (km/s)")


    # plot Vphi
    plt.subplot(2,2,2)
    plt.plot(seis_p/1.e9,mat_vphi_1/1.e3,color='r',linestyle='-',marker='^', \
             markerfacecolor='r',markersize=4)
    plt.plot(seis_p/1.e9,mat_vphi_2/1.e3,color='k',linestyle='-',marker='v', \
             markerfacecolor='k',markersize=4)
    plt.plot(seis_p/1.e9,mat_vphi_3/1.e3,color='b',linestyle='-',marker='x', \
             markerfacecolor='b',markersize=4)
    plt.plot(seis_p/1.e9,mat_vphi_4/1.e3,color='g',linestyle='-',marker='o', \
             markerfacecolor='g',markersize=4)
    plt.plot(seis_p/1.e9,mat_vphi_5/1.e3,color='y',linestyle='-',marker='*', \
             markerfacecolor='y',markersize=4)

    plt.title("Vphi (km/s)")

    # plot density
    plt.subplot(2,2,3)
    plt.plot(seis_p/1.e9,mat_rho_1/1.e3,color='r',linestyle='-',marker='^', \
             markerfacecolor='r',markersize=4,label='mgd')
    plt.plot(seis_p/1.e9,mat_rho_2/1.e3,color='k',linestyle='-',marker='v', \
             markerfacecolor='k',markersize=4,label='slb')
    plt.plot(seis_p/1.e9,mat_rho_3/1.e3,color='b',linestyle='-',marker='x', \
             markerfacecolor='b',markersize=4,label='slb3')
    plt.plot(seis_p/1.e9,mat_rho_4/1.e3,color='g',linestyle='-',marker='o', \
             markerfacecolor='g',markersize=4,label='bm2')
    plt.plot(seis_p/1.e9,mat_rho_5/1.e3,color='y',linestyle='-',marker='*', \
             markerfacecolor='y',markersize=4,label='bm3')
    plt.title("density (kg/m^3)")
    plt.legend(loc='upper left')


    plt.savefig("output_figures/example_compare_all_methods.png")
    plt.show()

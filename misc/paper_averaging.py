# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This script reproduces Cottaar, Heister, Rose, Unterborn (2014) Figure 2.

This example shows the effect of different averaging schemes. Currently four
averaging schemes are available:
1. Voight-Reuss-Hill
2. Voight averaging
3. Reuss averaging
4. Hashin-Shtrikman averaging

See Watt et al., 1976 Journal of Geophysics and Space Physics for explanations
of each averaging scheme.

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
import colors

if __name__ == "__main__":
    figsize=(6,5)
    prop={'size':12}
    #plt.rc('text', usetex=True)
    plt.rc('font', family='sanserif')
    figure=plt.figure(dpi=100,figsize=figsize)

    """ choose 'slb2' (finite-strain 2nd order shear modulus,
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

    amount_perovskite = 0.6
    method = 'slb3'

    rock = burnman.Composite( [ (minerals.SLB_2011.mg_perovskite(), amount_perovskite),
                    (minerals.SLB_2011.wuestite(), 1.0-amount_perovskite) ] )
    rock.set_method(method)

    perovskitite = burnman.Composite( [ (minerals.SLB_2011.mg_perovskite(), 1.0), ] )
    perovskitite.set_method(method)

    periclasite = burnman.Composite( [ (minerals.SLB_2011.wuestite(), 1.0), ] )
    periclasite.set_method(method)

    #seismic model for comparison:
    # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    seismic_model = burnman.seismic.PREM()
    #set on how many depth slices the computations should be done
    number_of_points = 20
    # we will do our computation and comparison at the following depth values:
    depths = np.linspace(700e3, 2800e3, number_of_points)
    #alternatively, we could use the values where prem is defined:
    #depths = seismic_model.internal_depth_list()
    pressures, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

    temperatures = burnman.geotherm.brown_shankland(pressures)


    print "Calculations are done for:"
    rock.debug_print()

        #calculate the seismic velocities of the rock using a whole battery of averaging schemes:

        # do the end members, here averaging scheme does not matter (though it defaults to Voigt-Reuss-Hill)
    rho_pv, vp_pv, vs_pv, vphi_pv, K_pv, G_pv = \
            burnman.velocities_from_rock(perovskitite, pressures, temperatures)
    rho_fp, vp_fp, vs_fp, vphi_fp, K_fp, G_fp = \
            burnman.velocities_from_rock(periclasite, pressures, temperatures)

        #Voigt Reuss Hill averaging
    rho_vrh, vp_vrh, vs_vrh, vphi_vrh, K_vrh, G_vrh = \
            burnman.velocities_from_rock(rock, pressures, temperatures, averaging_scheme=burnman.averaging_schemes.VoigtReussHill())

        #Voigt averaging
    rho_v, vp_v, vs_v, vphi_v, K_v, G_v = \
            burnman.velocities_from_rock(rock, pressures, temperatures, averaging_scheme=burnman.averaging_schemes.Voigt())

        #Reuss averaging
    rho_r, vp_r, vs_r, vphi_r, K_r, G_r = \
            burnman.velocities_from_rock(rock, pressures, temperatures, averaging_scheme=burnman.averaging_schemes.Reuss())

        #Upper bound for Hashin-Shtrikman averaging
    rho_hsu, vp_hsu, vs_hsu, vphi_hsu, K_hsu, G_hsu = \
            burnman.velocities_from_rock(rock, pressures, temperatures, averaging_scheme=burnman.averaging_schemes.HashinShtrikmanUpper())

        #Lower bound for Hashin-Shtrikman averaging
    rho_hsl, vp_hsl, vs_hsl, vphi_hsl, K_hsl, G_hsl = \
            burnman.velocities_from_rock(rock, pressures, temperatures, averaging_scheme=burnman.averaging_schemes.HashinShtrikmanLower())

    #linear fit
    vs_lin = vs_pv*amount_perovskite + vs_fp*(1.0-amount_perovskite)

    # PLOTTING

    # plot vs
    ax = figure.add_subplot(1,1,1)
    plt.plot(pressures/1.e9,vs_v/1.e3,color=colors.color(0),linewidth=2,linestyle='-',marker='^',\
        markersize=4,label='Voigt')
    plt.plot(pressures/1.e9,vs_r/1.e3,color=colors.color(5),linewidth=2,linestyle='-',marker='v',\
        markersize=4,label='Reuss')
    plt.plot(pressures/1.e9,vs_vrh/1.e3,color=colors.color(1),linestyle='-',marker='*',\
        markersize=6,label='Voigt-Reuss-Hill')
    plt.fill_between(pressures/1.e9, vs_hsu/1.e3, vs_hsl/1.e3, facecolor='red', lw=0, label='asdf',interpolate=False)




    #plt.plot(pressures/1.e9,vs_hsu/1.e3,color='r',linestyle='-',\
    #    markersize=4,label='Hashin-Shtrikman')
    #plt.plot(pressures/1.e9,vs_hsl/1.e3,color='r',linestyle='-',marker='x',\
    #    markersize=4)
    plt.plot(pressures/1.e9,vs_lin/1.e3,color='k',linewidth=2,linestyle='--',\
        markersize=4,label='linear')
    plt.plot(pressures/1.e9,vs_pv/1.e3,color=colors.color(2),linewidth=2,linestyle='-',marker='d',\
        markersize=4,label='Mg Perovskite')
    plt.plot(pressures/1.e9,vs_fp/1.e3,color=colors.color(4),linewidth=2,linestyle='-',marker='x',\
        markersize=6,label=r'W\"ustite')

    plt.ylim(3.0,7.5)
    plt.xlim(min(pressures)/1.e9,max(pressures)/1.e9)

    simArtist = plt.Line2D((0,1),(0,0), color='r', lw=5, linestyle='-')
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:3]+[simArtist]+handles[3:], labels[0:3]+['Hashin-Shtrikman']+labels[3:], loc='lower right',ncol=2,prop=prop)


    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Shear velocity $V_s$ (km/s)')
    plt.savefig("example_averaging.pdf",bbox_inches='tight')
    plt.show()

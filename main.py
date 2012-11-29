#system libs:
import numpy
import scipy.optimize as opt
import math
import matplotlib.pyplot as pyplot

#own libs:
import geotherm
import prem
from tools import *
from eos_from_ian import birch_murnaghan
import seismic

# TODO: add up weight percent and check <100 and tell them how much

molar_mass = {'Fe':55.845, 'Mg':24.305, 'O':15.999, 'Al':26.982, 'Ca':40.078, 'Si':28.085} # g/mol
Av = 6.022141e23 # Avogadro constant in 1/mol 
boltzmann_constant = 1.3806503e-23 # in m^2 kg s^-2 K^-1
gas_constant = Av * boltzmann_constant # in J mol^-1 K^-1

lower_mantle_mass = 4.043e27*.75 # in g




# convert weight percentage (amount, 1.00 = 100%) of a given element to molar mass
def weight_pct_to_mol(element, amount):

    return amount * lower_mantle_mass / molar_mass[element] * Av


def test_mol_conv():
    assert weight_pct_to_mol('Fe', 1.0) == 2*weight_pct_to_mol('Fe', 0.5)
    #assert float_eq(weight_pct_to_mol('Fe', 1.0), 3.26987875846e+49)


def conv_inputs(inp):
    names = {'Mg':'MgO','Fe':'FeO','Si':'SiO2', 'Ca':'Ca', 'Al':'Al'}
    out = {}
    for a in inp:
        out[names[a]] = weight_pct_to_mol(a,inp[a])
    return out
    

# compute phases of pv, fp, st
# inp = {'MgO':beta, 'FeO': , 'SiO2': gamma, 'Ca':, 'Al':} in mol
# params = {'Fe in pv': , 'Ca in pv':, 'Al in pv', 'Fe in fp':}
# returns: 'mol pv' A, 'mol fp' B, 'mol st' C in mol
# 'Mg in pv':0, 'Fe in pv':0, 'Ca in pv':0,'Si in pv':0, 'Al in pv':0
# 'Mg in fp':0,'Fe in fp':0
def determine_phases(inp, params):

    ret = {'mol pv':0., 'mol fp':0., 'mol st':0.}
    ret['Mg in pv'] = 1-params['Fe in pv']-params['Ca in pv'] 
    ret['Fe in pv'] = params['Fe in pv']
    ret['Ca in pv'] = params['Ca in pv']
    ret['Si in pv'] = 1-params['Al in pv']
    ret['Al in pv'] = params['Al in pv']
    ret['Mg in fp'] = 1 - params['Fe in fp']
    ret['Fe in fp'] = params['Fe in fp']
 
    beta = inp['MgO']
    gamma = inp['SiO2']

    if (beta > gamma):
        ret['mol pv'] = beta - gamma
        ret['mol fp'] = beta - ret['mol pv']
        ret['mol st'] = 0.
    elif (beta < gamma):
        ret['mol pv'] = beta
        ret['mol fp'] = 0.
        ret['mol st'] = gamma - ret['mol pv']
    else:
        ret['mol pv'] = beta
        ret['mol fp'] = 0.
        ret['mol st'] = 0.
    
    return ret



# test some composition (Javoy 2010, Table 6, PLoM)
inp1 = {'Mg':0.213, 'Fe': 0., 'Si':0.242, 'Ca':0., 'Al':0.} # wt%
inp2 = conv_inputs(inp1)
params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
t = determine_phases(inp2, params)
print inp1
print inp2
print t



def test_phases():
    # test everything into pv
    inp = {'MgO':20., 'FeO': 0., 'SiO2':20, 'CaO':0, 'Al2O3':0.}
    params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
    t = determine_phases(inp, params)
    assert t['mol pv'] == 20.
    assert t['mol fp'] == 0.
    assert t['mol st'] == 0.

    #
    inp = {'MgO':10, 'FeO': 0., 'SiO2':0., 'CaO':0., 'Al2O3':0.0}
    params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
    t = determine_phases(inp, params)
    assert t['mol pv'] == 10.
    assert t['mol fp'] == 0.
    assert t['mol st'] == 0.

    #
    inp = {'MgO':10, 'FeO': 0., 'SiO2':3, 'CaO':0., 'Al2O3':0.0}
    params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
    t = determine_phases(inp, params)
    assert t['mol pv'] == 7.
    assert t['mol fp'] == 3.
    assert t['mol st'] == 0.

    #
    inp = {'MgO':3., 'FeO': 0., 'SiO2':7., 'CaO':0., 'Al2O3':0.0}
    params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
    t = determine_phases(inp, params)
    assert t['mol pv'] == 3.
    assert t['mol fp'] == 0.
    assert t['mol st'] == 4.



#input: pv, fp, st in mol
#return: bulk modulus, shear modulus, density
def eqn_of_state(inp):
    # placeholder for now
    bla = 2.0

    out = {}
    out['density']= lambda pressure: 1+bla*pressure

    return out








#murakami test:

def compute_moduli(p,T,V0,K_0,K_prime,dKdT,a_0,a_1,gamma_0,molar_weight,atoms_per_unit_cell,n,eta_0s,q,G_0,G_prime):
    K_T = K_0 + dKdT*(T-300)
    alpha = a_0 + a_1*T
    P_th = alpha * K_T*(T-300)
    func = lambda x: birch_murnaghan (V0/x, 1., K_0, K_prime) + P_th - p
    V = opt.brentq(func, 0.1, 3.*V0)
    density = molar_weight*atoms_per_unit_cell / (Av*V*1e-24)    #correct according to jc and cayman

        # low T:
        # C_v = 234. * n * Av * boltzmann_constant (T/theta_D) ^ 3.  (3.43 from Poirier/EarthInterior)
        # high T:
        # C_v = 3. * n * gas_constant
        # n = number of moles
    C_v = 3. * n * gas_constant # in J * mol / K

    a2_s = -2.*gamma_0 - 2.*eta_0s # eq 47 (Stixrude)
    f = 1./2. * ( pow(V0/V ,2./3.) - 1.) # eq 24
    a1_ii = 6. * gamma_0 # eq 47
    a2_iikk = -12.*gamma_0+36.*pow(gamma_0,2.) - 18.*q*gamma_0 # eq 47
           
    nu_o_nu0_sq = 1.+ a1_ii*f + 1./2.*a2_iikk * pow(f,2.) # eq 41
    gamma = gamma_0 * pow(V0/V,-q)   # where from??
    eta_s = - gamma - 1./2. * pow( nu_o_nu0_sq, 2.) * pow(2.*f+1.,2.)*a2_s # eq 46
        
    #G = G_0 + G_prime*p + dGdT * (T-300.) #simple model
    
    delta_Uq = C_v* (T-300.) *1e3 * Av / lower_mantle_mass / 1e9  #cayman, others say 'what is this?'
    G = pow(1.+2.*f, 5./2.) * (G_0 + (3.*K_0*G_prime - 5.*G_0)*f \
                                   + (6.*K_0*G_prime - 24.*K_0 -14.*G_0 + 9./2.*K_0*K_prime)*pow(f,2.)) \
                                   - eta_s*density*delta_Uq #eq 33
    shear_mod = G

    # simple model:
    #K = K_0 + K_prime * p + dKdT*(T-300.)
    #K = K * (1. + alpha * gamma_0 * T)        # formula Matas, D6
    K = pow(1.+2.*f, 5./2.) * ( K_0 + (3*K_0*K_prime -5.*K_0)*f+27./2.*(K_0*K_prime-4.*K_0)*pow(f,2.)) \
        + (gamma+1.-q)*gamma * density * delta_Uq \
        - pow(gamma,2.) * density * C_v * (T-300.)*1e3 * Av / lower_mantle_mass / 1e9   # eq 32
    bulk_mod = K

    return V, density, bulk_mod, shear_mod


#input molar_abundance for pv and fp
#output: list_p, list_Vs, list_Vp, pv_density, fp_density, pv_shearmod, fp_shearmod, prem_shearmod
def murakami(molar_abundance):
    al_p = 0.075
    pv_X_Mg = 0.94
    fp_X_Mg = 0.79
    molar_weight=[pv_X_Mg*molar_mass['Mg']+(1.-pv_X_Mg)*molar_mass['Fe']+(1.-al_p)*molar_mass['Si']+al_p*molar_mass['Al']+3.*molar_mass['O'], \
                  fp_X_Mg*molar_mass['Mg']+(1.-fp_X_Mg)*molar_mass['Fe']+molar_mass['O']]

    list_p = []
    list_Vs = []
    list_Vp = []
    pv_density = []
    fp_density = []
    pv_shearmod = []
    fp_shearmod = []
    prem_shearmod = []

    for p in numpy.arange(28,141,5.):
        #T=geotherm.geotherm_formula(p)
        #T=geotherm.geotherm(p)
        T=geotherm.geotherm_brown(p) #by far the best fit

        density = [0., 0.]
        bulk_mod =[0., 0.]
        shear_mod = [0., 0.]

        # pv:
        # values from ricolleau table 1
        V0 = 164.
        K_0 = 245.
        K_prime = 4.
        dKdT = -.036
        a_0 = 3.19e-5
        a_1 = 0.88e-8

        gamma_0 = 1.48
        atoms_per_unit_cell = 4.
        G_0 = 166.
        G_prime = 1.57
        #dGdT = -.02 unused
        eta_0s = 2.4
        q = 1.4
        n = 5.

        #murakami supp. material Table 5:
        K_0 = 281.
        K_prime = 4.1
        G_0 = 173. #-25. #looks good :-)
        G_prime = 1.56

        V, density[0], bulk_mod[0], shear_mod[0] \
            = compute_moduli(p,T,V0,K_0,K_prime,dKdT,a_0,a_1,gamma_0, molar_weight[0], atoms_per_unit_cell,n, eta_0s, q, G_0, G_prime)
        pv_density.append(density[0])


        # fp:
        # values from ricolleau table 1
        if (p<=50.):
            V0 = 76.44
            K_0 = 158.
            K_prime = 4.
            dKdT = -.034
            a_0 = 2.20e-5
            a_1 = 3.61e-8
        else:
            V0 = 74.04
            K_0 = 170.
            K_prime = 4.
            dKdT = -.034
            a_0 = 2.20e-5
            a_1 = 3.61e-8
            #missing values, should we use 0?
            #dKdT = 0.
            #a_0 = 0.
            #a_1 = 0.

        gamma_0 = 1.50
        atoms_per_unit_cell = 4.
        if (p>=50.): # reading from fig 3 in Murakami (X_Mg=0.79)
            G_0 = 116.  # low spin
            G_prime = 1.65
            #dGdT = -.02 unused
        else:
            G_0 = 103.  # high spin
            G_prime = 1.78
            #dGdT = -.02 unused
        eta_0s = 3.0
        q = 1.5
        n = 2.

        V, density[1], bulk_mod[1],shear_mod[1] \
            = compute_moduli(p,T,V0,K_0,K_prime,dKdT,a_0,a_1,gamma_0, molar_weight[1], atoms_per_unit_cell,n, eta_0s, q, G_0, G_prime)
        fp_density.append(density[1])



        #if (p>=50.): # from the text in Murakami (page 2), X_Mg = 0.92
        #    shear_mod[1] = 130. + 2.04*p -.02*(T-300) # low spin
        #else:
        #    shear_mod[1] = 113. + 2.15*p -.02*(T-300)



        pv_shearmod.append(shear_mod[0])
        fp_shearmod.append(shear_mod[1])

        V_p,V_s,V_phi = seismic.get_velocities(molar_abundance, molar_weight, bulk_mod, shear_mod, density, T)
        list_p.append(p)
        list_Vp.append(V_p)
        list_Vs.append(V_s)

        # shearmodulus = density * V_s^2: (Sanne)
        prem_shearmod.append(prem.prem_density(p)*pow(prem.prem_V(p)[1],2.0))

    return list_p, list_Vs, list_Vp, pv_density, fp_density, pv_shearmod, fp_shearmod, prem_shearmod


#compute prem
prem_p = numpy.arange(28.3,135.0,5)
prem_vp = [prem.prem_V(y)[0] for y in prem_p]
prem_vs = [prem.prem_V(y)[1] for y in prem_p]
prem_density = [prem.prem_density(y) for y in prem_p]

#compute murakami for 100% fp
molar_abundance=[0., 1.0]
list_p, fp_Vs, fp_Vp, pv_density, fp_density, pv_shearmod, fp_shearmod, prem_shearmod \
    = murakami(molar_abundance)

molar_abundance=[1.0, .0]
_, pv_Vs, pv_Vp, _,_,_,_,_ \
    = murakami(molar_abundance)

#molar_abundance=[0.95, .05]
molar_abundance=[0.93, .07]
_, mix_Vs, mix_Vp, _,_,_,_,_ \
    = murakami(molar_abundance)

mix_density = [molar_abundance[0] * pv_density[i] + molar_abundance[1] * fp_density[i] for i in range(len(pv_density))]

# plot Vs
pyplot.subplot(2,2,1)
p1,=pyplot.plot(list_p,fp_Vs,'-k')
p2,=pyplot.plot(list_p,pv_Vs,'-b')
p3,=pyplot.plot(list_p,mix_Vs,'-r')
p4,=pyplot.plot(prem_p,prem_vs,'--k',markerfacecolor='white')
pyplot.legend([p1,p2,p3,p4],["fp", "pv", "mix (pv: "+str(molar_abundance[0]*100.)+"%)", "PREM"], loc=4)
pyplot.title("Vs")
pyplot.xlim(25,135)
pyplot.ylim(5.,7.6)

# plot Vp
pyplot.subplot(2,2,2)
p1,=pyplot.plot(list_p,fp_Vp,'-k')
p2,=pyplot.plot(list_p,pv_Vp,'-b')
p3,=pyplot.plot(list_p,mix_Vp,'-r')
p4,=pyplot.plot(prem_p,prem_vp, '--k',markerfacecolor='white')
pyplot.legend([p1,p2,p3],["fp", "pv", "mix", "PREM"], loc=4)
pyplot.title("Vp")
pyplot.xlim(30,135)
pyplot.ylim(9.25,14.)

# plot shear mod
#pyplot.subplot(2,2,4)
#pyplot.title("Shearmodulus comparison")
#p1,=pyplot.plot(list_p,fp_shearmod,'-g')
#p2,=pyplot.plot(list_p,pv_shearmod,'-b')
#p3,=pyplot.plot(list_p,prem_shearmod,'--k',markerfacecolor='white',markevery=1)
#pyplot.legend([p1,p2,p3],["fp", "pv", "PREM"], loc=4)
#pyplot.xlim(30,135)

# plot density
pyplot.subplot(2,2,3)
p1,=pyplot.plot(list_p,fp_density,'-k')
p2,=pyplot.plot(list_p,pv_density,'-b')
p3,=pyplot.plot(prem_p,prem_density,'--k',markerfacecolor='white')
p4,=pyplot.plot(list_p,mix_density,'-r')
pyplot.legend([p1,p2,p3,p4],["fp", "pv", "PREM", "mix"], loc=4)
pyplot.title("density")
pyplot.xlim(30,135)
pyplot.ylim(4.,6.5)



pyplot.show()




test_phases()
test_mol_conv()












#print "full example:"

#inp1 = {'Mg':0.5, 'Fe': 0, 'Si':0.5, 'Ca':0.0, 'Al':0} # wt%
#inp2 = conv_inputs(inp1)
#print "in:", inp1
#print "out:", inp2

#params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
#t = determine_phases(inp2, params)
#print "phases:", t

#ret = eqn_of_state(t)
#
#print "eos:", ret
#
#print ret['density'](42)









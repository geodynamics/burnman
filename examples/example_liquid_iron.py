import os, sys, numpy as np, matplotlib.pyplot as plt
import matplotlib.image as mpimg
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))
import burnman

from scipy.optimize import fsolve

atomic_masses=burnman.processchemistry.read_masses()

class liquid_iron( burnman.Mineral ):
    def __init__(self):
        formula='Fe1.0'
        formula = burnman.processchemistry.dictionarize_formula(formula)
        m = burnman.processchemistry.formula_mass(formula, atomic_masses)
        rho_0 = 7019.
        V_0 = m/rho_0
        D = 7766.
        Lambda = 1146.
        self.params = {
            'name': 'liquid iron',
            'formula': formula,
            'equation_of_state': 'aa',
            'P_0': 1.e5,
            'T_0': 1811.,
            'S_0': 100., # to fit
            'molar_mass': m,
            'V_0': V_0,
            'E_0': 0.,
            'K_S': 109.7e9,
            'Kprime_S': 4.661,
            'Kprime_prime_S': -0.043e-9,
            'grueneisen_0': 1.735,
            'grueneisen_prime': -0.130/m*1.e-6,
            'grueneisen_n': -1.870,
            'a': [248.92*m, 289.48*m],
            'b': [0.4057*m, -1.1499*m],
            'Theta': [1747.3, 1.537],
            'theta': 5000.,
            'lmda': [-325.23*m, 302.07*m, 30.45*m],
            'xi_0': 282.67*m,
            'F': [D/rho_0, Lambda/rho_0],
            'n': sum(formula.values()),
            'molar_mass': m}
        burnman.Mineral.__init__(self)


        
liq = liquid_iron()

burnman.tools.check_eos_consistency(liq, P=10.e9, T=7000., tol=0.01, verbose=True)

# Gibbs
P0 = 10.e9
T0 = 3000.
dP = 10.
dT = 1.
liq.set_state(P0, T0)
G0 = liq.gibbs
S0 = liq.S
V0 = liq.V
H0 = liq.H
F0 = liq.helmholtz
E0 = liq.internal_energy
Cv0 = liq.heat_capacity_v


P1 = liq.method.pressure(T0+dT , liq.V, liq.params)
liq.set_state(P1, T0 + dT)

# Check Cv, S, E, V calculations
Cv1 = liq.heat_capacity_v
S1 = liq.S
E1 = liq.internal_energy

print 'Cv', Cv0, (E1 - E0)/dT, T0*(S1 - S0)/dT
print 'dEdP|V', (E1 - E0)/(P1 - P0), liq.V/liq.gr
print 'dEdS|V', (E1 - E0)/(S1 - S0), T0 + 0.5*dT

def diffS(args, T1, P0, T0, S0, m):
    P1 = args[0]
    liq.set_state(P0, T0)
    S0 = liq.S
    liq.set_state(P1, T1)
    S1 = liq.S
    return S1 - S0


def diffV(args, T1, P0, T0, V0, m):
    P1 = args[0]
    liq.set_state(P0, T0)
    V0 = liq.V
    liq.set_state(P1, T1)
    V1 = liq.V
    return V1 - V0

T1 = T0 + 1.
P1 = fsolve(diffS, [P0], args=(T1, P0, T0, S0, liq))[0]
E1 = liq.internal_energy
V1 = liq.V
print 'dEdV|S', (E1 - E0)/(V0 - V1)/1.e9, P0/1.e9


# Check isentropic energy change
dV = V0*1.e-5
E0i = liq.method._isentropic_energy_change(V0, liq.params)
E1i = liq.method._isentropic_energy_change(V0+dV, liq.params)
Pi = liq.method._isentropic_pressure(V0, liq.params)
print 'Pth', -(E1i - E0i)/dV/1.e9, Pi/1.e9



# P constant
P1 = P0
T1 = T0 + dT
liq.set_state(P0, T0+dT)
G1 = liq.gibbs
F1 = liq.helmholtz
H1 = liq.H
S1 = liq.S

# T constant
P2 = P0 + dP
T2 = T0
liq.set_state(P0+dP, T0)
G2 = liq.gibbs
H2 = liq.H
V2 = liq.V
F2 = liq.helmholtz

# S constant
T3 = T0+dT
P3 = fsolve(diffS, [P0], args=(T0 + dT, P0, T0, V0, liq), xtol=1.e-20)[0]
liq.set_state(P3, T0+dT)
H3 = liq.H
E3 = liq.internal_energy
V3 = liq.V

# V constant
T4 = T0 + dT
P4 = fsolve(diffV, [P0], args=(T1, P0, T0, V0, liq), xtol=1.e-20)[0]
F4 = liq.helmholtz
E4 = liq.internal_energy
S4 = liq.S

print 'dG = -SdT + VdP'
print -(G1 - G0)/dT, (S1 + S0)/2. # P constant
print (G2 - G0)/dP, (V2 + V0)/2. # T constant

print 'dH = TdS + VdP'
print (H1 - H0)/(S1 - S0), (T1 + T0)/2. # P constant
print (H3 - H0)/(P3 - P0), (V2 + V0)/2. # S constant

print 'dE = TdS - PdV'
print (E4 - E0)/(S4 - S0), (T4 + T0)/2. # V constant
print -(E3 - E0)/(V3 - V0)/1.e9, (P3 + P0)/2./1.e9  # S constant

print 'dF = -SdT - PdV'
print -(F4 - F0)/(T4 - T0), (S4 + S0)/2. # V constant
print -(F2 - F0)/(V2 - V0)/1.e9, (P2 + P0)/2./1.e9 # T constant


exit()


# Find heat capacities
temperatures = np.linspace(1000., 15000., 101)
Cvs = np.empty_like(temperatures)
m = 0.055845
rhos = np.empty_like(temperatures)
densities = [5.e3,10.e3, 15.e3]
#densities = [8.5e3, 9.56e3, 10.81e3, 12.28e3, 14.03e3, 16.14e3, 18.68e3]
for rho in densities:
    V = m/rho
    for i, T in enumerate(temperatures):
        Cvs[i] = liq.method.heat_capacity_v(0., T, V, liq.params)/burnman.constants.gas_constant
        #Cvs[i] = liq.method._C_v_el(V, T, liq.params)/burnman.constants.gas_constant

    plt.plot(temperatures, Cvs)

fig1 = mpimg.imread('../burnman/data/input_figures/AA1994_liq_iron_TCv_different_densities.png')
plt.imshow(fig1, extent=[1000., 15000., 0., 6.], aspect='auto')
plt.ylim(0., 6.)
plt.title('AA1994, Figure 5')
plt.show()




def temperature(T, P, rho, mineral):
    mineral.set_state(P, T[0])
    return mineral.density - rho

# Check grueneisen values
Prhos = [[1.e5, 7019.],
         [0.2e9, 5500.],
         [0.2e9, 6000.],
         [0.2e9, 6500.],
         [277.4e9, 12643.],
         [331.5e9, 13015.],
         [397.1e9, 13417.]]
for Prho in Prhos:
    P, rho = Prho
    T = fsolve(temperature, 1811., args=(P, rho, liq))[0]

    K0 = liq.K_T
    dP = 100.
    liq.set_state(P+dP, T)
    dK = liq.K_T - K0

    grueneisen_Irvine_Stacey_1975 =  (0.5*dK/dP - 5./6. + 2./9.*P/K0)/(1. - 4.*P/(3.*K0))
    print P, rho, T, liq.grueneisen_parameter, grueneisen_Irvine_Stacey_1975


liq.set_state(130.e9, 5000.)
print liq.gr

# Find volumes and temperatures up the reference isentrope
# Check standard state values
liq.set_state(1.e5, 1811.)
print liq.gr, liq.alpha, liq.K_S, liq.C_p, liq.density
reference_entropy = liq.S

def isentrope(T, P, Sref, mineral):
    mineral.set_state(P, T[0])
    return Sref - mineral.S

pressures = np.linspace(0., 500.e9, 101)
temperatures = np.empty_like(pressures)
rhos = np.empty_like(pressures)
Vps = np.empty_like(pressures)

for i, P in enumerate(pressures):
    temperatures[i] = fsolve(isentrope, 1811., args=(P, reference_entropy, liq))[0]
    rhos[i] = liq.density
    Vps[i] = np.sqrt(liq.K_S/liq.density)

np.savetxt(header='Pressures (GPa)   Temperatures (K)   Densities (kg/m^3), Vps (km/s)', X=zip(*[pressures/1.e9, temperatures, rhos, Vps]), fname='isentropic_PTrhoVp.dat')

fig1 = mpimg.imread('../burnman/data/input_figures/AA1994_liq_iron_PTrho_reference_isentrope.png')

plt.imshow(fig1, extent=[0.0, 500., 6., 15.], aspect='auto')
plt.plot(pressures/1.e9, rhos/1.e3, linewidth=2) 
plt.title('AA1994 Figure B1 (1/2)')
plt.show()

plt.imshow(fig1, extent=[0.0, 500., 1500., 7000.], aspect='auto')
plt.plot(pressures/1.e9, temperatures, linewidth=2)
plt.title('AA1994 Figure B1 (2/2)')
plt.show()


# Find densities at 1 bar
temperatures = np.linspace(1800., 2400., 100)
rhos = np.empty_like(temperatures)
rhos_mizuno = np.empty_like(temperatures)
Vps = np.empty_like(temperatures)


P = 1.e5
for i, T in enumerate(temperatures):
    liq.set_state(1.e5, T)
    rhos[i] = liq.density/1.e3
    Vps[i] = np.sqrt(liq.p_wave_velocity)
    rhos_mizuno[i] = (7162. - 0.735*(T - 1808.))/1.e3

plt.plot(temperatures, Vps)
plt.title('Vp as a function of T at 1 bar')
plt.show()

fig1 = mpimg.imread('../burnman/data/input_figures/AA1994_liq_iron_Trho_1bar.png')
plt.imshow(fig1, extent=[1800., 2400., 6.65, 7.1], aspect='auto')
plt.plot(temperatures, rhos)
plt.plot(temperatures, rhos_mizuno)
plt.ylim(6.65, 7.1)
plt.title('AA1994 Figure 1')
plt.show()


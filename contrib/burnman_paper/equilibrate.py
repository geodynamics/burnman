import matplotlib.pyplot as plt
import numpy as np

import burnman_path  # adds the local burnman directory to the path
from burnman import minerals
import burnman
assert burnman_path  # silence pyflakes warning
from burnman import equilibrate

bdg = minerals.SLB_2011.mg_fe_bridgmanite()
ppv = minerals.SLB_2011.post_perovskite()
per = minerals.SLB_2011.ferropericlase()
cpv = minerals.SLB_2011.ca_perovskite()
cf = minerals.SLB_2011.ca_ferrite_structured_phase()

bdg.set_composition([0.85, 0.1, 0.05])
ppv.set_composition([0.85, 0.1, 0.05])
per.set_composition([0.7, 0.3])
cf.set_composition([0.4, 0.3, 0.3])
assemblage = burnman.Composite([bdg, per, cpv, cf])

KLB1 = burnman.Composition({'SiO2': 44.48,
                            'Al2O3': 3.59,
                            'FeO': 8.10,
                            'MgO': 39.22,
                            'CaO': 3.44,
                            'Na2O': 0.30}, 'weight')


KLB1.renormalize('atomic', 'total', 1.)
composition = KLB1.atomic_composition
print(composition)

equality_constraints = [('P', 30.e9), ('T', 2000.)]

sol, prm = equilibrate(composition, assemblage, equality_constraints)
assert sol.success
print(sol.assemblage)

print(sol.assemblage.S*sol.assemblage.n_moles, sol.assemblage.n_moles)

print(bdg.molar_fractions)
print(per.molar_fractions)
#bdg.set_composition([0.85, 0.1, 0.05])
per.set_composition([0.7, 0.3])
#cf.set_composition([0.4, 0.3, 0.3])


assemblage = burnman.Composite([bdg, per, cpv, cf, ppv])
assemblage.set_state(100.e9, 2000.)
equality_constraints = [('S', sol.assemblage.S*sol.assemblage.n_moles),
                        ('phase_fraction', (ppv, 0.))]


sol, prm = equilibrate(composition, assemblage, equality_constraints)
assert sol.success
print(sol.assemblage)

assemblage = burnman.Composite([bdg, per, cpv, cf, ppv])
assemblage.set_state(120.e9, 2600.)
equality_constraints = [('S', sol.assemblage.S*sol.assemblage.n_moles),
                        ('phase_fraction', (bdg, 0.))]


sol, prm = equilibrate(composition, assemblage, equality_constraints)
assert sol.success
print(sol.assemblage)
print(np.linalg.cond(sol.J))


J = sol.J.copy()

print(prm.parameter_names)
fx = np.ones(14)
fF = np.ones(14)
for i, name in enumerate(prm.parameter_names):
    if name == 'Pressure (Pa)':
        fx[i] = 1.e9
    elif name == 'Temperature (K)':
        fx[i] = 1.e2
    elif name[:2] == 'x(':
        fx[i] = 0.05
    elif name[:2] == 'p(':
        fx[i] = 0.05

fF[0] = 1. # entropy
fF[1] = 0.05 # phase fraction
fF[2:8] = 1000. # reactions
fF[8:] = 0.001 # bulk
J2 = np.einsum("ij, i, j ->ij", J, 1./fF, fx)
#print(J2)
print(f'{np.linalg.cond(J2):.2e}')

"""

equality_constraints = [['S', 53.4366987322636],
                        ['X', [np.array([0.,  0.,  1.,  0.,
                                      0., -0.,  0., -0., -0.,
                                      0.,  0., -0.,  0., 0.]), 0.0]]]

F = burnman.tools.equilibration.F(sol.x, sol.assemblage, equality_constraints,
                                  prm.reduced_composition_vector)

m = 1.e-8
dx = np.eye(14)*m
J_test = np.zeros((14, 14))
for i in range(14):
    Fi = burnman.tools.equilibration.F(sol.x + dx[i], sol.assemblage, equality_constraints,
                                       prm.reduced_composition_vector)
    J_test[:,i] = (Fi - F)/m
assert sol.success
print(sol.assemblage)

np.set_printoptions(precision=3)
print(J_test - sol.J)
print(np.max(J_test - sol.J))
imax = np.argmax(J_test - sol.J)
print(imax%14, int((imax - imax%14)/14))
print((J_test - sol.J)[int((imax - imax%14)/14),imax%14])
print(prm.parameter_names[imax%14])
"""

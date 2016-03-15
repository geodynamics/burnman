# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import
import numpy as np
from .material import Material
from .averaging_schemes import AveragingScheme


class Model(object):

    """
    planetary model realization

    initialize with a rock, averaging scheme and list of pressures and temperatures

    then access any of the members to get lists of properties for pressures and temperatures

    all computations are done automatically and lazily

    """

    def __init__(self, rock, p, T, avgscheme):
        assert(len(p) == len(T))
        assert(len(p) > 0)
        assert(isinstance(rock, Material))
        assert(isinstance(avgscheme, AveragingScheme))

        self.rock = rock
        self.p = p
        self.T = T
        self.avgscheme = avgscheme

        self.moduli = None
        self.mat_rho = None
        self.mat_vs = None
        self.mat_vp = None
        self.mat_vphi = None
        self.mat_K = None
        self.mat_G = None
        self.mat_vp = None
        self.mat_vs = None
        self.mat_vphi = None
        self.alpha = None
        self.c_p = None
        self.c_v = None
        self.mat_V = None

    def v_s(self):
        self.compute_velocities_()
        return self.mat_vs

    def v_p(self):
        self.compute_velocities_()
        return self.mat_vp

    def v_phi(self):
        self.compute_velocities_()
        return self.mat_vphi

    def density(self):
        self.avg_moduli_()
        return self.mat_rho

    def K(self):
        self.avg_moduli_()
        return self.mat_K

    def G(self):
        self.avg_moduli_()
        return self.mat_G

    def thermal_expansivity(self):
        if self.alpha is None:
            self.calc_moduli_()
            n_pressures = len(self.p)
            self.alpha = np.zeros(n_pressures)
            for idx in range(n_pressures):
                V = np.array([m['V'] for m in self.moduli[idx]])
                alpha = np.array([m['alpha'] for m in self.moduli[idx]])
                self.alpha[
                    idx] = self.avgscheme.average_thermal_expansivity(V, alpha)

        return self.alpha

    def heat_capacity_p(self):
        self.calc_heat_capacities_()
        return self.c_p

    def heat_capacity_v(self):
        self.calc_heat_capacities_()
        return self.c_v

    def calc_moduli_(self):
        """
        Internal function to compute the moduli if necessary.
        """
        if self.moduli is None:
                self.moduli = [[] for p in self.p]

                for idx in range(len(self.p)):
                    self.rock.set_state(self.p[idx], self.T[idx])
                    (minerals, fractions) = self.rock.unroll()
                    for (mineral, fraction) in zip(minerals, fractions):
                        e = {}
                        e['fraction'] = fraction
                        e['V'] = fraction * mineral.molar_volume
                        e['K'] = mineral.adiabatic_bulk_modulus
                        e['G'] = mineral.shear_modulus
                        e['rho'] = mineral.molar_mass / mineral.molar_volume
                        e['alpha'] = mineral.thermal_expansivity
                        e['c_v'] = mineral.heat_capacity_v
                        e['c_p'] = mineral.heat_capacity_p
                        self.moduli[idx].append(e)

    def avg_moduli_(self):
        """
        Internal function to average moduli if necessary.
        """
        if self.mat_V is None:
            self.calc_moduli_()
            n_pressures = len(self.p)
            self.mat_V = np.zeros(len(self.p))
            self.mat_K = np.zeros(len(self.p))
            self.mat_G = np.zeros(len(self.p))
            self.mat_rho = np.zeros(len(self.p))

            for idx in range(n_pressures):
                mods = self.moduli[idx]

                fractions = np.array([e['fraction'] for e in mods])
                V_frac = np.array([m['V'] for m in mods])
                K_ph = np.array([m['K'] for m in mods])
                G_ph = np.array([m['G'] for m in mods])
                rho_ph = np.array([m['rho'] for m in mods])

                self.mat_V[idx] = sum(V_frac)
                self.mat_K[idx] = self.avgscheme.average_bulk_moduli(
                    V_frac, K_ph, G_ph)
                self.mat_G[idx] = self.avgscheme.average_shear_moduli(
                    V_frac, K_ph, G_ph)
                self.mat_rho[idx] = self.avgscheme.average_density(
                    V_frac, rho_ph)

    def calc_heat_capacities_(self):
        """
        Internal function to compute the heat capacities if necessary.
        """
        if self.c_p is None:
            self.calc_moduli_()
            self.c_v = np.zeros(len(self.p))
            self.c_p = np.zeros(len(self.p))
            for idx in range(len(self.p)):
                fractions = np.array([m['fraction'] for m in self.moduli[idx]])
                alphas = np.array([m['alpha'] for m in self.moduli[idx]])
                c_v = np.array([m['c_v'] for m in self.moduli[idx]])
                c_p = np.array([m['c_p'] for m in self.moduli[idx]])
                self.c_v[idx] = self.avgscheme.average_heat_capacity_v(
                    fractions, c_v)
                self.c_p[idx] = self.avgscheme.average_heat_capacity_p(
                    fractions, c_p)

    def compute_velocities_(self):
        """
        Internal function to compute the velocities if necessary.
        """
        if self.mat_vp is None:
            self.avg_moduli_()
            self.mat_vs = np.ndarray(len(self.p))
            self.mat_vp = np.ndarray(len(self.p))
            self.mat_vphi = np.ndarray(len(self.p))

            for i in range(len(self.p)):
                self.mat_vs[i] = np.sqrt(self.mat_G[i] / self.mat_rho[i])
                self.mat_vp[i] = np.sqrt(
                    (self.mat_K[i] + 4. / 3. * self.mat_G[i]) / self.mat_rho[i])
                self.mat_vphi[i] = np.sqrt(self.mat_K[i] / self.mat_rho[i])

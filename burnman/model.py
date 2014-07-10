# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
from burnman.main import compute_velocities
import main


class Model:
    """
    planetary model realization

    initialize with a rock, averaging scheme and list of pressures and temperatures

    then access any of the members to get lists of properties for pressures and temperatures

    all computations are done automatically and lazily

    """
    def __init__(self, rock, p, T, avgscheme):
        self.rock = rock
        self.p = p
        self.T = T
        self.avgscheme = avgscheme

        self.moduli = None
        self.avgmoduli = None
        self.mat_rho = None
        self.mat_vs = None
        self.mat_vp = None
        self.mat_vphi = None
        self.mat_K = None
        self.mat_G = None
        self.mat_vp = None
        self.mat_vs = None
        self.mat_vphi = None

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
        if self.mat_rho is None:
            self.avg_moduli_()
            self.mat_rho = np.array([m.rho for m in self.avgmoduli])
        return self.mat_rho

    def K(self):
        if self.mat_K is None:
            self.avg_moduli_()
            self.mat_K = np.array([m.K for m in self.avgmoduli])
        return self.mat_K

    def G(self):
        if self.mat_G is None:
            self.avg_moduli_()
            self.mat_G = np.array([m.G for m in self.avgmoduli])
        return self.mat_G


    def calc_moduli_(self):
        if self.moduli is None:
            self.moduli = main.calculate_moduli(self.rock, self.p, self.T)

    def avg_moduli_(self):
        if self.avgmoduli is None:
            self.calc_moduli_()
            self.avgmoduli = main.average_moduli(self.moduli, self.avgscheme)

    def compute_velocities_(self):
        if self.mat_vp is None:
            self.avg_moduli_()
            self.mat_vp, self.mat_vs, self.mat_vphi = main.compute_velocities(self.avgmoduli)








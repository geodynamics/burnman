import unittest
import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from util import BurnManTest
from burnman.nonlinear_fitting import *


class test_fitting(BurnManTest):

    def test_linear_fit(self):
        
        i, x, Wx, y, Wy = np.loadtxt('../burnman/data/input_fitting/Pearson_York.dat', unpack=True)

        data = np.array([x, y]).T
        cov = np.array([[1./Wx, 0.*Wx], [0.*Wy, 1./Wy]]).T

        class m():
            def __init__(self, data, cov, guessed_params, delta_params):
                self.data = data
                self.data_covariances = cov
                self.set_params(guessed_params)
                self.delta_params = delta_params
                self.mle_tolerances = np.array([1.e-1] * len(data[:,0])) # irrelevant for a linear model

            def set_params(self, param_values):
                self.params = param_values
        
            def get_params(self):
                return self.params
    
            def function(self, x, flag):
                return np.array([x[0], self.params[0]*x[0] + self.params[1]])

            def normal(self, x, flag):
                n = np.array([self.params[0], -1.])
                return n/np.linalg.norm(n)

        guessed_params = np.array([-0.5, 5.5])
        delta_params = np.array([1.e-3, 1.e-3]) # unimportant for a linear model
        fitted_curve = m(data, cov, guessed_params, delta_params)
        nonlinear_least_squares_fit(model=fitted_curve,
                                    param_tolerance = 1.e-5)
        
        self.assertArraysAlmostEqual([fitted_curve.WSS], [11.8663531941])

    def test_fit_PVT_data(self):
        fo = burnman.minerals.HP_2011_ds62.fo()

        pressures = np.linspace(1.e9, 2.e9, 8)
        temperatures = np.ones_like(pressures) * fo.params['T_0']

        PTV = np.empty((len(pressures), 3))

        for i in range(len(pressures)):
            fo.set_state(pressures[i], temperatures[i])
            PTV[i] = [pressures[i], temperatures[i], fo.V]

        params = ['V_0', 'K_0', 'Kprime_0']
        fitted_eos = burnman.eos_fitting.fit_PTV_data(fo, params, PTV, verbose=False)
        zeros = np.zeros_like(fitted_eos.pcov[0])

        self.assertArraysAlmostEqual(fitted_eos.pcov[0], zeros)

if __name__ == '__main__':
    unittest.main()

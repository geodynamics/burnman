from __future__ import absolute_import
import unittest
import os
import sys

sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals
from burnman import seismic

from util import BurnManTest


class test_seismic(BurnManTest):

    def test_internal_depth_list(self):
        models = [burnman.seismic.PREM(), burnman.seismic.STW105(),
                  burnman.seismic.AK135(), burnman.seismic.IASP91()]

        ref_depth_lists = {'PREM': [0.0, 6371000.0, 94],
                           'STW105': [0.0, 6371000.0, 750],
                           'AK135': [0.0, 6371000.0, 145],
                           'IASP91': [0.0, 6371000.0, 152]}

        for model in models:
            dl = model.internal_depth_list()
            name = model.__class__.__name__
            stats = [min(dl), max(dl), len(dl)]
            # print model.__class__.__name__, stats
            self.assertArraysAlmostEqual(stats, ref_depth_lists[name])

    def test_evaluate(self):
        models = [burnman.seismic.PREM(),
                  burnman.seismic.Fast(),
                  burnman.seismic.Slow(),
                  burnman.seismic.STW105(), burnman.seismic.AK135(),
                  burnman.seismic.IASP91()]

        ref = {
            'PREM': [12817.6924, 6932.8549000000003, 10010.358087588364, 5120.6290999999992],
            'Fast': [12795.360742611414, 6941.4225568201909, 9973.8053779446909, 5120.6290999999992],
            'Slow': [12795.360742611414, 6904.3291880013148, 10008.07546360332, 5120.6290999999992],
            'STW105': [12817.85012987013, 6932.9342065251822, 10010.486818122406, 5120.6749255622426],
            'AK135': [12798.468686868688, 6920.5212121212116, 9997.1320353488773, 5103.8969696969698],
            'IASP91': [12794.4, 6921.0, 9991.4805389391604], }

        for model in models:
            dl = model.internal_depth_list()
            name = model.__class__.__name__
            depth = 2000e3
            vars = ['v_p', 'v_s', 'v_phi', 'density']  # skip gravity
            if name == "IASP91":
                vars = vars[0:-1]  # skip density
            result = model.evaluate(vars, [depth])
            result = list(result.T[0])
            # print "'%s': %s," % (name, result)
            self.assertArraysAlmostEqual(result, ref[name])


if __name__ == '__main__':
    unittest.main()

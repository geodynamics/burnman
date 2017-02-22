from __future__ import absolute_import
from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

import argparse
import os.path
import sys
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))
    
from burnman.perplex import create_perplex_table

parser = argparse.ArgumentParser(description='Call werami to create a burnman-readable tab file.')

parser.add_argument('--werami_path', metavar='path', type=str, nargs='+', required=True,
                    help='The path to werami')
parser.add_argument('--project', metavar='project name', type=str, nargs='+', required=True,
                    help='The name of the project file (without the suffix)')
parser.add_argument('--outfile', metavar='output file', type=str, nargs='+', required=True,
                    help='The name of the output table file')
parser.add_argument('--n_pressures', type=int, nargs=1, required=True,
                    help='The number of pressure steps in the grid')
parser.add_argument('--n_temperatures', type=int, nargs=1, required=True,
                    help='The number of pressure steps in the grid')
parser.add_argument('--pressure_range', type=float, nargs=2,
                    help='Minimum and maximum values of pressure (Pa; optional)')
parser.add_argument('--temperature_range', type=float, nargs=2,
                    help='Minimum and maximum values of temperature (K; optional)')

args = parser.parse_args()
if not hasattr(args, 'pressure_range'):
    args.pressure_range = None
if not hasattr(args, 'temperature_range'):
    args.temperature_range = None

create_perplex_table(args.werami_path, args.project[0], args.outfile[0], args.n_pressures[0], args.n_temperatures[0], args.pressure_range, args.temperature_range)


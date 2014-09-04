# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

# This is a standalone program that converts a tabulated version of the Stixrude and Lithgow-Bertelloni data format into the standard burnman format (printed to stdout)


import sys


def read_dataset(datafile):
    f=open(datafile,'r')
    ds=[]
    for line in f:
        ds.append(line.decode('utf-8').split())
    return ds

ds=read_dataset('data/raw_endmember_datasets/slb_2011.txt')

print 'from burnman.mineral import Mineral'
print 'from burnman.processchemistry import read_masses, dictionarize_formula, formula_mass'
print ''
print 'atomic_masses=read_masses(\'data/input_masses/atomic_masses.dat\')'
print ''

formula='0'
for idx, m in enumerate(ds):
    if idx == 0:
        param_names=m
    else:
        print 'class', m[0].lower(), '(Mineral):'
        print '    def __init__(self):'
        print ''.join(['       formula=\'',m[1],'\''])
        print '       formula = dictionarize_formula(formula)'
        print '       self.params = {'
        print ''.join(['            \'name\': \'', m[0], '\','])
        print '            \'formula\': formula,'
        print '            \'equation_of_state\': \'slb3\','
        for pid, param in enumerate(m):
            if pid > 1 and pid%2 == 0:
                print '            \''+param_names[pid]+'\':', float(param), ','
        print '            \'n\': sum(formula.values()),'
        print '            \'molar_mass\': formula_mass(formula, atomic_masses)}'
        print ''
        print '       self.uncertainties = {'
        for pid, param in enumerate(m):
            if pid > 1 and pid%2 == 1 and pid<21:
                print '            \''+param_names[pid]+'\':', float(param), ','
        pid=21
        param=m[pid]
        print '            \''+param_names[pid]+'\':', float(param), '}'
        print ''

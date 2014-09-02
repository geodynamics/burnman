# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

# This is a standalone program that converts the Holland and Powell data format into the standard burnman format (printed to stdout)
# It only outputs properties of solid endmembers - other endmembers are currently ignored.


import sys

# Components
components=['Si','Ti','Al','Fe','Mg','Mn','Ca','Na','K','O','H','C','Cl','e-','Ni','Zr','S','Cu','Cr']
class Endmember:
    def __init__(self,name,atoms,formula,sites,comp,H,S,V,Cp,a,k,flag,od):
        self.name=name # Name of end member
        self.atoms=atoms # Number of atoms in the unit formula
        self.formula=formula
        self.sites=sites # Notional number of sites
        self.comp=comp # Composition
        self.H=H # Enthalpy
        self.S=S # Entropy
        self.V=V # Volume
        self.Cp=Cp # Heat capacity (c0 + c1*T + c2*T**2 + c3*T**(-0.5)) 
        self.a=a # Thermal expansion
        self.k=k # Bulk Modulus (and first two derivatives wrt pressure)
        self.flag=flag
        self.od=od

def read_dataset(datafile):
    f=open(datafile,'r')
    ds=[]
    for line in f:
        ds.append(line.split())
    return ds

ds=read_dataset('data/raw_endmember_datasets/tc-ds62.txt')

def getmbr(ds, mbr):
    mbrarray = []
    for i in range(0, int(ds[0][0])):
        if ds[i*4+3][0] == mbr:
            atoms = 0.0
            formula=''
            for j in range(3,len(ds[i*4+3])-1,2):
                atoms = atoms + float(ds[i*4+3][j])
                formula=formula+components[int(ds[i*4+3][j-1])-1]+str(round(float(ds[i*4+3][j]),10))
            if mbr.endswith('L'):
                flag=-2
                od=[0]
            else:
                flag=int(ds[i*4+6][4])
            endmember=Endmember(mbr,atoms,formula, int(ds[i*4+3][1]), map(float,ds[i*4+3][2:(len(ds[i*4+3])-1)]), float(ds[i*4+4][0]), float(ds[i*4+4][1]), float(ds[i*4+4][2]), map(float,ds[i*4+5]), float(ds[i*4+6][0]), map(float,ds[i*4+6][1:4]), flag, map(float,ds[i*4+6][5:]))
            return endmember

print 'from burnman.mineral import Mineral'
print 'from burnman.processchemistry import read_masses, dictionarize_formula, formula_mass'
print ''
print 'atomic_masses=read_masses(\'data/input_masses/atomic_masses.dat\')'
print ''

formula='0'
for i in range(int(ds[0][0])):
    mbr=ds[i*4+3][0]
    M=getmbr(ds,mbr)
    if mbr == 'and': # change silly abbreviation
        mbr = 'andalusite'
    if M.flag != -1 and M.flag != -2 and M.k[0] > 0:
        print 'class', mbr, '(Mineral):'
        print '    def __init__(self):'
        print ''.join(['       formula=\'',M.formula,'\''])
        print '       formula = dictionarize_formula(formula)'
        print '       self.params = {'
        print ''.join(['            \'name\': \'', M.name, '\','])
        print '            \'formula\': formula,'
        print '            \'equation_of_state\': \'mtait\','
        print '            \'H_0\':', M.H*1e3, ','
        print '            \'S_0\':', M.S*1e3, ','
        print '            \'V_0\':', M.V*1e-5, ','
        print '            \'Cp\':', [round(M.Cp[0]*1e3,10), round(M.Cp[1]*1e3,10), round(M.Cp[2]*1e3,10), round(M.Cp[3]*1e3,10)], ','
        print '            \'a_0\':', M.a, ','
        print '            \'K_0\':', M.k[0]*1e8, ','
        print '            \'Kprime_0\':', M.k[1], ','
        print '            \'Kdprime_0\':', M.k[2]*1e-8, ','
        print '            \'n\': sum(formula.values()),'
        if M.flag==0:
            print '            \'molar_mass\': formula_mass(formula, atomic_masses)}'
        else:
            print '            \'molar_mass\': formula_mass(formula, atomic_masses),'
        if M.flag==1:
            print '            \'landau_Tc\':', M.od[0], ','
            print '            \'landau_Smax\':', M.od[1]*1e3, ','
            print '            \'landau_Vmax\':', M.od[2]*1e-5, '}'
        if M.flag==2:
            print '            \'BW_deltaH\':', M.od[0]*1e3, ','
            print '            \'BW_deltaV\':', M.od[1]*1e-5, ','
            print '            \'BW_W\':', M.od[2]*1e3, ','
            print '            \'BW_Wv\':', M.od[3]*1e-5, ','
            print '            \'BW_n\':', M.od[4], ','
            print '            \'BW_factor\':', M.od[5], '}'

        print ''
print ''

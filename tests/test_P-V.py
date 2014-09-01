import os, sys
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from   burnman import minerals

P0=2.e10
T=2273.15

wad=minerals.other.Katsura_2009_wadsleyite()
wad.set_method("slb3")

# Here we find the volume at the P, T of interest
wad.set_state(P0,T)
print 'V @', T, 'K,', P0, 'Pa:'
print wad.V, 'm^3/mol'
print ''

# Here we take the V from the previous step and calculate the pressure. Hopefully it's the same as our original input!!
V=wad.V

P1=wad.eos_pressure(T,V)
print 'Retrieved pressure:', P1, 'Pa'
print 'Pressure difference:', P1-P0, 'Pa'
print 'Fractional difference:', (P1-P0)/P0

"""
    BurnMan- a lower mantle toolkit
    Copyright (C) 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""


import numpy as np
import matplotlib.pyplot as plt

K=251
Kp=4.1
G=175
Gp=np.linspace(1.5,2.3, 50)

saveratio=np.empty_like(Gp)

for i in range(len(Gp)):
  term1=3.*K*Gp[i]-5.*G
  term2=6.*K*Gp[i]-24.*K-14.*G+9./2.*K*Kp
  saveratio[i]=term2*0.01/(term1*0.1)
  print Gp[i], term1, term2, term2/term1, term2*0.01/(term1*0.1)

term3=3.*K*Gp[25]-3.*K-7.*G
term4=-9.*K*Kp-3.*G*Gp[25]+30.*K+63.*G
print term4/term3


plt.plot(Gp,saveratio)

plt.show()

# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

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

# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2019 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
example_holland_powell_uncertainties
------------------------------------

This extremely short example script shows how one can visualize
and manipulate uncertainties in zero-point energies as found in the
Holland and Powell dataset.

*Uses:*

* :doc:`mineral_database`


*Demonstrates:*

* creating basic layer
* calculating thermoelastic properties with self-consistent pressures
* seismic comparison
"""
from __future__ import absolute_import

# Here we import standard python modules that are required for
# usage of BurnMan.
import numpy as np
import matplotlib.pyplot as plt


import burnman_path  # adds the local burnman directory to the path
# Here we import the relevant modules from BurnMan.
from burnman.minerals import HP_2011_ds62
from burnman.optimize.nonlinear_fitting import plot_cov_ellipse  # for plotting

assert burnman_path  # silence pyflakes warning

plt.style.use('ggplot')

# First, we read in the covariance matrix
cov = HP_2011_ds62.cov()

# Now we can find the desired rows and columns of the covariance matrix
# by cross-referencing against the list of mineral names
# (quartz, periclase, enstatite)
indices = [cov['endmember_names'].index(name) for name in ['q', 'per', 'en']]

# Slice the required rows and columns from the covariance matrix
Cov_selected = cov['covariance_matrix'][np.ix_(indices,indices)]

# The following line transforms the covariance matrix so that we can look
# at the uncertainties associated with the endmember reaction
# quartz + periclase = 0.5*enstatite
A = np.array([[1., 1., 0.],
              [0., 0., 0.5]])
cov_transformed = A.dot(Cov_selected).dot(A.T)

sigma_x = np.sqrt(cov_transformed[0][0])
sigma_y = np.sqrt(cov_transformed[1][1])
corr_xy = cov_transformed[0][1]/sigma_x/sigma_y

print('sigma(q+per) = {0:.2f} J/mol'.format(sigma_x))
print('sigma(en/2) = {0:.2f} J/mol'.format(sigma_y))
print('corr(q+per,en/2) = {0:.2f}'.format(corr_xy))


# Finally, we plot the covariance matrix
fig = plt.figure(figsize = (5, 5))
ax = fig.add_subplot(1, 1, 1)

plot_cov_ellipse(cov_transformed, [0., 0.], nstd=1., ax=ax)

ax.set_xlim(-500., 500.)
ax.set_ylim(-500., 500.)
ax.set_xlabel('$\sigma_{q+per}$ (J/mol MgSiO$_3$)')
ax.set_ylabel('$\sigma_{en/2}$ (J/mol MgSiO$_3$)')
plt.show()

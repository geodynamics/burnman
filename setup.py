from __future__ import absolute_import
import re

versionstuff = dict(
    re.findall("(.+) = '(.+)'\n", open('burnman/version.py').read()))

metadata = dict(name='burnman',
                version=versionstuff['version'],
                description=('A thermoelastic and thermodynamic toolkit '
                             'for Earth and planetary sciences'),
                url='https://geodynamics.github.io/burnman/',
                author='The BurnMan Team',
                author_email='bob.myhill@bristol.ac.uk',
                license='GPL',
                long_description_content_type="text/x-rst",
                long_description=('BurnMan is a Python library for generating '
                                  'thermodynamic and thermoelastic models of '
                                  'materials and planetary interiors.'),
                packages=['burnman',
                          'burnman.utils',
                          'burnman.minerals',
                          'burnman.eos',
                          'burnman.calibrants',
                          'burnman.tools',
                          'burnman.classes',
                          'burnman.optimize'],
                package_data={'burnman': ['data/input_*/*']},
                classifiers=['License :: OSI Approved :: GNU General Public '
                             'License v2 or later (GPLv2+)',
                             'Programming Language :: Python :: 3.7',
                             'Programming Language :: Python :: 3.8',
                             'Programming Language :: Python :: 3.9',
                             'Programming Language :: Python :: 3.10'],
                )

# Try to use setuptools in order to check dependencies.
# if the system does not have setuptools, fall back on
# distutils.
try:
    from setuptools import setup
    metadata['install_requires'] = ['numpy', 'matplotlib', 'scipy', 'sympy']
except ImportError:
    from distutils.core import setup


setup(**metadata)

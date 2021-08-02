from __future__ import absolute_import
import re

versionstuff = dict(
    re.findall("(.+) = '(.+)'\n", open('burnman/version.py').read()))

metadata = dict(name='burnman',
                version=versionstuff['version'],
                description='a thermoelastic and thermodynamic toolkit for Earth and planetary sciences',
                url='http://burnman.org',
                author='The BurnMan Team',
                author_email='bob.myhill@bristol.ac.uk',
                license='GPL',
                long_description='BurnMan is a Python library for generating thermodynamic and thermoelastic models of planetary interiors.',
                packages=['burnman', 'burnman.minerals', 'burnman.eos'],
                package_data={'burnman': ['data/input_*/*']},
                classifiers=[
                'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
                'Programming Language :: Python :: 3.6',
                'Programming Language :: Python :: 3.7',
                'Programming Language :: Python :: 3.8'],
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

from __future__ import absolute_import
import re

versionstuff = dict(
    re.findall("(.+) = '(.+)'\n", open('burnman/version.py').read()))

metadata = dict(name='burnman',
                version=versionstuff['version'],
                description='a thermoelastic and thermodynamic toolkit for Earth and planetary sciences',
                url='http://burnman.org',
                author='Ian Rose',
                author_email='ian.rose@berkeley.edu',
                license='GPL',
                long_description='BurnMan is a Python library for generating thermodynamic and thermoelastic models of planetary interiors.',
                packages=['burnman', 'burnman.minerals', 'burnman.eos'],
                package_data={'burnman': ['data/input_*/*']},
                classifiers=[
                'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
                'Programming Language :: Python :: 2.7',
                'Programming Language :: Python :: 3.4'],
                )

# Try to use setuptools in order to check dependencies.
# if the system does not have setuptools, fall back on
# distutils.
try:
    from setuptools import setup, find_packages
    metadata['install_requires'] = ['numpy', 'matplotlib', 'scipy']
except ImportError:
    from distutils.core import setup


setup(**metadata)

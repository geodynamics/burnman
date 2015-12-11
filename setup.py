from __future__ import absolute_import
import re

versionstuff = dict(re.findall("(.+) = '(.+)'\n", open('burnman/version.py').read()))

metadata = dict( name= 'burnman',
                 version = versionstuff['version'],
                 description='a lower mantle toolkit',
                 url='http://burnman.org',
                 author='Ian Rose',
                 author_email='ian.rose@berkeley.edu',
                 license='GPL',
                 long_description='BurnMan calculates elastic properties of lower mantle assemblages constrained by mineral physics.',
                 packages = ['burnman', 'burnman.minerals', 'burnman.eos'],
                 package_data = { 'burnman' : ['data/input_*/*'] },
               )

#Try to use setuptools in order to check dependencies.
#if the system does not have setuptools, fall back on
#distutils.
try:
  from setuptools import setup, find_packages
  metadata['install_requires'] = ['numpy', 'matplotlib', 'scipy']
except ImportError:
  from distutils.core import setup


setup ( **metadata )

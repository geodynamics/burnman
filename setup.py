from distutils.core import setup

setup ( name = 'burnman',
        version = '0.6',
        packages = ['burnman', 'burnman.minerals'],
        package_data = { 'burnman' : ['data/input_*/*'] }
       )

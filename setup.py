from setuptools import setup, find_packages

setup ( name = 'burnman',
        version = '0.6b3',
        description='a lower mantle toolkit',
        url='http://burnman.org',
        author='Ian Rose',
        author_email='ian.rose@berkeley.edu',
        license='GPL',
        long_description='BurnMan is a lower mantle shear velocity generator constrained by mineral physics.',
        packages = ['burnman', 'burnman.minerals'],
        package_data = { 'burnman' : ['data/input_*/*'] },
        install_requires=['numpy', 'matplotlib', 'scipy']
       )

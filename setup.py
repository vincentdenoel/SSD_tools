from setuptools import setup

setup(
    name='SSDpy',
    version='0.1',
    description='A collection of tools for the Structural and Stochastic Dynamics of systems',
    author='SSD research group at the University of Liège',
    author_email='ssd@uliege.be',
    #packages=['SSDpy','SSDpy.dyn','SSDpy.dyn.integrators','SSDpy.wind'],
    packages=['PackageName','PackageName.SubModule'],
    install_requires=['numpy', 'scipy', 'matplotlib'], # external packages as dependencies
    )
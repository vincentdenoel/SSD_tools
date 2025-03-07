from setuptools import setup, find_packages


setup(
    name='SSDpy',
    version='0.1',
    description='A collection of tools from the Structural and Stochastic Dynamics',
    author='SSD research group at the University of Liège - supervised by V. Denoël',
    author_email='ssd@uliege.be',
    packages=find_packages(),
    install_requires=['numpy', 'scipy', 'matplotlib'], # external packages as dependencies
    )


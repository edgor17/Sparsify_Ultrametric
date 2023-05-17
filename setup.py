from setuptools import setup

setup(
    name='Sparsify_Ultrametric',
    url='https://github.com/edgor17/Sparsify-Ultrametric',
    author='Evan Gorman',
    author_email='evan.gorman@colorado.edu',
    packages=['Sparsify_Ultrametric'],
    install_requires=['numpy','pandas','ete3','scipy','matplotlib'],
    version='0.1',
    description='Sparsification of Strictly Ultrametric Matrices',
    long_description=open('README.txt').read(),
)

from setuptools import setup, find_packages
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
setup(
    name='Sparsify_Ultrametric',
    url='https://github.com/edgor17/Sparsify-Ultrametric',
    author='Evan Gorman',
    author_email='evan.gorman@colorado.edu',
    packages=['Sparsify_Ultrametric'],
    install_requires=['numpy','pandas','ete3','scipy','matplotlib','scikit-learn'],
    version='0.1',
    description='Sparsification of Strictly Ultrametric Matrices',
    long_description=long_description
)

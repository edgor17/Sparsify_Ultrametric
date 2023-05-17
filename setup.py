from setuptools import setup
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
setup(
    name='SparsifyUltrametric',
    url='https://github.com/edgor17/Sparsify-Ultrametric',
    author='Evan Gorman',
    author_email='evan.gorman@colorado.edu',
    packages=['SparsifyUltrametric'],
    install_requires=['numpy','pandas','ete3','scipy','matplotlib'],
    version='0.1',
    description='Sparsification of Strictly Ultrametric Matrices',
    long_description=long_description
)

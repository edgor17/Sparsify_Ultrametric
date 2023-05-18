# Sparsify-Ultrametric

This package contains all data and scripts to reproduce the analysis in "Sparsification of Large Ultrametric Matrices: Insights into the Microbial Tree of Life"

To install:

'''
git clone https://github.com/edgor17/Sparsify_Ultrametric
cd Sparsify_Ultrametric/
pip install .
'''

The main goal is to take a binary tree structure and to produce the associated Haar-like wavelet basis and sparsified covariance matrix. We use ete3 to read .nwk formatted trees into a python data structure.  

'''
from Sparsify_Ultrametric import Sparsify
from ete3 import Tree

tree = Tree("/raw_data/97_otus_unannotated.tree",format=1)
[haarlike,pseudodiag]=Sparsify(tree)
'''

This computation may be slow, we have included the precomputed haarlike and pseudodiag files for convenience.


Given data consisting of OTU abundances it is necessary to map these abundances onto a reference phylogenetic tree. We require a table of OTU counts where each row is a sample and the columns are OTU ID's . The first row must consist of the OTU ids in quotations, these OTU ids are mapped to the leaves of the given tree.

example. py shows how to compute the Haar-like distances (and resulting PCoA embedding) of microbial mat samples from the Guerrero Negro hypersaline microbial mat (Harris et al., 2012). We also show how to obtain the spectrogram-like plots of Haar_like projections comparing two environments. 

Acknowledgments

The fast method for matching indices in SPARSIFY is thanks to a StackOverflow answer by Olivier Melançon found here: https://stackoverflow.com/questions/49247506/how-to-efficiently-find-the-indices-of-matching-elements-in-two-lists


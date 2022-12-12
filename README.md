# Sparsify-Ultrametric

Requires: Numpy, SciPy, ETE Tooklkit, sci-kit bio

Sparsify takes a newick tree structure as an input and produces its assocaited Haar-like wavelet basis and sparsified covariance matrix. Note that most phylogenetic trees (such as Greengenes) are actually two ORB-Trees connected to an external root with branch length zero. This is implicity assumed in Sparsify.

Given data consisting of OTU abundances it is necessary to map these abundances onto a reference phylogenetic tree. We require a table of OTU counts where each row is a sample and the columns are OTU ID's . The first row must consist of the OTU ids in quotations, these OTU ids are mapped to the leaves of the given tree.

example. py shows how to compute the Haar-like distances (and resulting PCoA embedding) of microbial mat samples from the Guerrero Negro hypersaline microbial mat (Harris et al., 2012). We also show how to obtain the spectrogram-like plots of Haar_like projections comparing two environments. 

Acknowledgments

The fast method for matching indices in SPARSIFY is thanks to a StackOverflow answer by Olivier Melançon found here: https://stackoverflow.com/questions/49247506/how-to-efficiently-find-the-indices-of-matching-elements-in-two-lists


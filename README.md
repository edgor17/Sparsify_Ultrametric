# Sparsify-Ultrametric

Requires: Numpy, SciPy, ETE Tooklkit, sci-kit bio

SPARSIFY takes a newick tree structure as an input and produces its assocaited Haar-like wavelet basis and sparsified covariance matrix. Note that most phylogenetic trees (such as Greengenes) are actually two ORB-Trees connected to an external root with branch length zero. This is implicity assumed in SPARSIFY.

haarlikeproj then projects OTU data over the Haar-like basis. This requires a table of OTU samples each row is an environment and each column the associated OTU counts. The first column must consist of the OTU ids in quotations, these OTU ids are mapped to the leaves of the given tree.

haarlikedistances computes the pairwise Haar-like distances between all environments. The resulting distance matrix can be exported for visualization in QIIME2 or a PCoA embedding can be computed using skbio. 


Acknowledgments
The fast method for matching indices in SPARSIFY is thanks to a StackOverflow answer by Olivier Melançon found here: https://stackoverflow.com/questions/49247506/how-to-efficiently-find-the-indices-of-matching-elements-in-two-lists


# Sparsify-Ultrametric

Requires: Numpy, SciPy, ETE Tooklkit, sci-kit bio

SPARSIFY takes a newick tree structure as an input and produces its assocaited Haar-like wavelet basis and sparsified covariance matrix. Note that most phylogenetic trees (such as Greengenes) are actually two ORB-Trees connected to an external root with branch length zero. This is implicity assumed in SPARSIFY.

haarlikeproj then projects OTU data over the Haar-like basis. This requires a table of OTU samples each row is an environment and each column the associated OTU counts. The first column must consist of the OTU ids in quotations, these OTU ids are mapped to the leaves of the given tree.

haarlikedistances computes the pairwise Haar-like distances between all environments. The resulting distance matrix can be exported for visualization in QIIME2 or a PCoA embedding can be computed using skbio. 



As an example we can compute the Haar-like distances (and resulting PCoA embedding) of microbial mat samples from the Guerrero Negro hypersaline microbial mat (Harris et al., 2012). We obtained this data through QIITA, which uses Greengenes 97% as their default phylogenetic tree to compare samples. The first step in our analysis is to compute the sparsified covariance matrix associated with Greengenes 97%. The original newick file is 97_otus_unannotated.tree and the Haar-like basis is 97haarlike.npz. The associated sparse covariance matrix is over 25mb and cannot be saved in this repository, however all that is required for the Haar-like distance is the diagonal of this matrix which is saved as eigest97.npy. 

Next the OTU data (microbialmatall.csv) is projected over the Haar-like basis producing a vector of scalars associated with each environment. These vectors are then scaled by eigest97 and summed to compute the Haar-like distances. Embedding these results using PCoA results in the following plot where the darker greys represent shallower soil samples:

![alt text](https://https://github.com/edgor17/Sparsify-Ultrametric/blob/main/soilgradient.png?raw=true "Title")

Acknowledgments

The fast method for matching indices in SPARSIFY is thanks to a StackOverflow answer by Olivier Melançon found here: https://stackoverflow.com/questions/49247506/how-to-efficiently-find-the-indices-of-matching-elements-in-two-lists


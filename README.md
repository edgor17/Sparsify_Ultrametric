# Sparsify-Ultrametric

This package contains all data and scripts to reproduce the analysis in "Sparsification of Large Ultrametric Matrices: Insights into the Microbial Tree of Life"

### Installation 

```
git clone https://github.com/edgor17/Sparsify_Ultrametric
cd Sparsify_Ultrametric/
pip install .
```

### Computation of the Haar-like Basis and Sparsified Covariance

The main goal is to take a binary tree structure and to produce the associated Haar-like wavelet basis and sparsified covariance matrix. We use ete3 to read .nwk formatted trees into a python data structure.  

```
from Sparsify_Ultrametric import Sparsify
from ete3 import Tree

tree = Tree("/raw_data/97_otus_unannotated.tree",format=1)
[haarlike,pseudodiag]=Sparsify(tree)
```

This computation may be slow, we have included the precomputed haarlike and pseudodiag files for convenience.

### The Haar-like Distance

To demonstrate the Haar-like distance we consider microbial mat samples from the Guerrero Negro hypersaline microbial mat (Harris et al., 2012). This data was obtained from [QIITA](https://qiita.ucsd.edu/study/description/1200#). We read in the feature table (OTU counts) and metadata (containing sample depths) then project these OTU counts onto the Haar-like basis resulting in "mags"

```
featuretable=pd.read_csv("/raw_data/mat/otus.txt", sep='\t')
metadata=pd.read_csv("/raw_data/metadata.txt", sep='\t')
[X,Y,mags,dic]=PreProcess(featuretable,metadata,'end_depth','regression',tree,haarlike)
```

To compute the Haar-like distances between these samples we need to rescale by lambda_v, obtained from the diagonal of our sparsified covariance matrix. 

```
lambdav=scipy.sparse.csr_matrix.diagonal(pseudodiag)
[D,modmags]=compute_Haar_dist(mags,lambdav)
```


Acknowledgments

The fast method for matching indices in SPARSIFY is thanks to a StackOverflow answer by Olivier Melançon found here: https://stackoverflow.com/questions/49247506/how-to-efficiently-find-the-indices-of-matching-elements-in-two-lists


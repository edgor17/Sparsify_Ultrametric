# Sparsify-Ultrametric

This package contains all data and scripts to reproduce the analysis in "Sparsification of Large Ultrametric Matrices: Insights into the Microbial Tree of Life"

For any questions email evan.gorman@colorado.edu

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

tree = Tree("Sparsify_Ultrametric/raw_data/97_otus_unannotated.tree",format=1)
[haarlike,pseudodiag]=Sparsify.sparsify(tree)
```

The sparsify computation may take a while to complete, we have included the precomputed haarlike and pseudodiag files for convenience, which can be loaded as:

```
import scipy
haarlike=scipy.sparse.load_npz('Sparsify_Ultrametric/precomputed/97haarlike.npz')
pseudodiag=scipy.sparse.load_npz('Sparsify_Ultrametric/precomputed/97pseudodiag.npz')
```

Additionally, all plots in the [paper](https://www.biorxiv.org/content/10.1101/2022.08.21.504697v2) can be reproduced using these precomputed files by following plots.py

### The Haar-like Distance

To demonstrate the Haar-like distance we consider microbial mat samples from the Guerrero Negro hypersaline microbial mat (Harris et al., 2012). This data was obtained from [QIITA](https://qiita.ucsd.edu/study/description/1200#). We read in the feature table (OTU counts) and metadata (containing sample depths) then project these OTU counts onto the Haar-like basis collecting the resulting magnitudes in "mags"

```
from Sparsify_Ultrametric.preprocessing import PreProcess
import pandas as pd

featuretable=pd.read_csv("Sparsify_Ultrametric/raw_data/otus.txt", sep='\t')
metadata=pd.read_csv("Sparsify_Ultrametric/raw_data/metadata.txt", sep='\t')
[X,Y,mags,dic]=PreProcess(featuretable,metadata,'end_depth','regression',tree,haarlike)
```

To compute the Haar-like distances between these samples we need to rescale by lambda_v, obtained from the diagonal of our sparsified covariance matrix. 

```
import scipy
from Sparsify_Ultrametric.utils import compute_Haar_dist

lambdav=scipy.sparse.csr_matrix.diagonal(pseudodiag)
[D,modmags]=compute_Haar_dist(mags,lambdav)
```

We can then plot the associated PCoA embedding.

```
from Sparsify_Ultrametric.utils import PCoA
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np

depths=Y.values.astype(float)
haarlikecoord=PCoA(D,2)
color = iter(cm.viridis(np.linspace(0, 1, 9)))
fig,ax=plt.subplots(1,figsize=(8, 3))
for item in np.unique(depths):
    c=next(color)
    haarlikeplot=ax.scatter(-haarlikecoord[np.where(depths==item),0] ,haarlikecoord[np.where(depths==item),1],c=c,label=str(item/1000))
ax.title.set_text('Haar-like Distance Embedding')
plt.legend(title="Depth in Meters",loc='center left', bbox_to_anchor=(1.1, .5))
```

![haarlikeex](https://github.com/edgor17/Sparsify_Ultrametric/assets/87628022/3304d8a8-fc6a-4195-b6a4-657e44075a9e)




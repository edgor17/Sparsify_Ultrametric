import numpy as np
import scipy as scipy
from sklearn.manifold import MDS

#The computation of the sparsified matrix takes a few hours (for large trees). However, this only needs to be done once ever for
#any tree of interest. To save time in this example we've precomputed both the Haar-like basis as "97haarlike.npz" and the diagonal 
#of the sparsified matrix as "eigest.npy". To compute this from scratch use pseudodiag=sparsify("treepath") followed by eigest=pseudodiag.diagonal()

eigest=np.load("eigest97.npy")
Haar_like=scipy.sparse.load_npz("97haarlike.npz")

#load the OTU abundance data from the microbial mat study
data=np.transpose(np.loadtxt("/Users/Evan/Desktop/pythontrees/microbialmatall.csv",delimiter=','))

#map the OTU counts onto the 97 Greengenes reference phylogeny
abundvecs=Match_to_tree(data,"97_otus_unannotated.tree")

#compute the pairwise Haar-like distances between environments. This is done by projecting OTU abundance onto the Haar-like 
#vectors and rescaling by the estimated eigenvalues of the 97 Greengenes covariance matrix
[D,mags]=compute_Haar_dist(abundvecs,Haar_like,eigest,normalized=True)

#convert mags to dense format for plotting
mags=mags.todense()

#Spectrogram-like plot comparing the first and second environments w.r.t the Haar-like basis. The sum of the magnitude of the 
#spikes in this plot is the Haar-like distance between the two environments. Accordingly, this plot shows how each split in the 
#phylogeny contributed to the resulting Haar-like distance. 
plt.plot(np.square(mags[0,:]-mags[1,:]).T)

#Compute the coordinates of the PCoA embedding of D
embedding = MDS(n_components=2,dissimilarity='precomputed')
X_transformed = embedding.fit_transform(D)

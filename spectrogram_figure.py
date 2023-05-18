import numpy as np
import scipy
from ete3 import Tree
import matplotlib.pyplot as plt
from Sparsify_Ultrametric import utils

#load in precomputed lambda_v and Haar-like basis 
eigest=np.load("eigest97.npy")
Haar_like=scipy.sparse.load_npz("97haarlike.npz")

#load a table of OTU counts
data=np.transpose(np.loadtxt("/Users/Evan/Desktop/pythontrees/microbialmatall.csv",delimiter=','))

#load the .nwk tree
tree = Tree("97_otus_unannotated.tree",format=1)

#map the OTU counts onto the reference tree
abundvecs=utils.Match_to_tree(data,tree)

#Project the OTU abundances onto the Haar-like basis and compute the Haar-like distance between all pairs of samples
[D,mags]=utils.compute_Haar_dist(abundvecs,Haar_like,eigest,normalized=True)

#Convert the resulting projections to dense format for plotting purposes
mags=mags.todense()

#Plot the spectrogram comparing the 
plt.plot(np.square(mags[0,:]-mags[1,:]).T)

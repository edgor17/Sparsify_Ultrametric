import pandas as pd
import scipy as scipy
import scipy.sparse
import numpy as np
from ete3 import Tree

def compute_Haar_dist(data,Haar_like,weightvec,normalized):
    #data should consist of OTU abundance counts where each row is a sample and each column is an OTU ID. 
    #The weightvec weights projections onto each Haar-like wavelet. Using the diagonal of the sparsified phylogenetic covariance 
    #produces distances similar to that of DPCOA. normalized=True will divide each row by the sum of its OTU counts. 
    #The output is a pairwise distance matrix that can be visualized via PCoA or some other embedding method. 
    abunds=scipy.sparse.csr_matrix(data)
    if normalized:
        abunds=abunds/abunds.sum(0)
        abunds=scipy.sparse.csr_matrix(abunds)
    mags=Haar_like@abunds
    modmags=np.transpose(np.asarray(mags.todense())* weightvec[:, np.newaxis])
    
    N=len(data[0,:])
    
    #Build Haar-like distance matrix
    D=np.zeros((N,N))
    for i in range(N):
        print(i)
        for j in range(i+1,N):
            distdiff=((modmags[i,:]-modmags[j,:]))
            d=np.sum(np.square(distdiff))
            D[i,j]=np.sqrt(d)        
    D=D+np.transpose(D)
    return D

def Match_to_tree(data, tree):
    #maps the given OTU abundances onto a reference tree. Data must include OTU ids as columns 
    #tree should be a filepath to a .nwk file
    t = Tree(tree,format=1)
    allleaves=t.get_leaves()
    leafcount=len(t)
    abundvec=np.zeros((data.shape[0]-1,leafcount))
    N=data.shape[0]-1



    #Index leaves of the tree to find OTUs
    leaflist=[]
    for i in range(leafcount):
        leaf=allleaves[i]
        leaflist.append(leaf.name)

    #match OTUs to leaves in the tree    
    for i in range(data.shape[1]):
        num=data[0,i].astype(int)
        name=num.astype(str)
        leafloc=leaflist.index(name)
        for j in range(1,N+1):
            abundvec[j-1,leafloc]=data[j,i]

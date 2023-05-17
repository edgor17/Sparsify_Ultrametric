import pandas as pd
import scipy as scipy
import numpy as np
from ete3 import Tree

def compute_Haar_dist(abundvecs,Haar_like,weightvec,normalized):
    '''

    Parameters
    ----------
    data : np.ndarray
        OTU abundance counts where each row is a sample and each column is an OTU ID
    Haar_like : scipy.sparse.csr_matrix
        Haar_like basis vectors
    weightvec : np.ndarray
        lamba_v's to use in distance computation
    normalized : bool
        True if abundance counts should be normalized to sum to 1

    Returns
    -------
    D : np.ndarray
        Pairwise Haar-like distances
    modmags : scipy.sparse.csr_matrix
        rescaled (by weightvec) projections of OTU abundance onto Haar-like basis

    '''
    
    abunds=scipy.sparse.csr_matrix(abundvecs.T)
    if normalized:
        abunds=abunds/abunds.sum(0)
        abunds=scipy.sparse.csr_matrix(abunds)
    mags=Haar_like@abunds
    modmags=np.transpose(np.asarray(mags.todense())* np.sqrt(weightvec[:, np.newaxis]))
    modmags=scipy.sparse.csr_matrix(modmags)
    
    N=len(abundvecs[:,0])
    
    #Build Haar-like distance matrix
    D=np.zeros((N,N))
    for i in range(N):
        for j in range(i+1,N):
            distdiff=((modmags[i,:]-modmags[j,:]))
            d=scipy.sparse.csr_matrix.sum(scipy.sparse.csr_matrix.power(distdiff,2))
            D[i,j]=np.sqrt(d)        
    D=D+np.transpose(D)
    return D, modmags

def Match_to_tree(data, tree):
    '''
    Parameters
    ----------
    data : np.ndarray
        OTU abundanes. Must include OTU ids as columns
    tree : TreeNode 
        Tree

    Returns
    -------
    abundvec : np.ndarray
        OTU abundances as mapped onto the tree

    '''
    
    allleaves=tree.get_leaves()
    leafcount=len(tree)
    abundvecs=np.zeros((data.shape[0]-1,leafcount))
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
            abundvecs[j-1,leafloc]=data[j,i]
    return abundvecs

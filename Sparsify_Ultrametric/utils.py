import pandas as pd
import scipy as scipy
import numpy as np
from ete3 import Tree

def compute_Haar_dist(mags,weightvec):
    '''

    Parameters
    ----------
    mags : np.ndarray
        Projection of OTU abundance count onto Haar-like basis
    weightvec : np.ndarray
        lamba_v's to use in distance computation

    Returns
    -------
    D : np.ndarray
        Pairwise Haar-like distances
    modmags : scipy.sparse.csr_matrix
        rescaled (by weightvec) projections of OTU abundance onto Haar-like basis

    '''
    
    modmags=np.transpose(np.asarray(mags.todense())* np.sqrt(weightvec[:, np.newaxis]))
    modmags=scipy.sparse.csr_matrix(modmags)
    
    N=mags.shape[1]
    
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

def PCoA(D,n):
    A2=np.square(D)
    J=np.eye(len(D))-1/(len(D))*np.ones(len(D))
    B=-1/2*J@A2@J
    [D,V]=np.linalg.eigh(B)
    D=np.real(D)
    D=np.flip(D)
    V=np.flip(V,axis=1)
    V=V[:,0:n]
    D=np.sqrt(D[0:n])
    X=V@np.diag(D)
    
    return X

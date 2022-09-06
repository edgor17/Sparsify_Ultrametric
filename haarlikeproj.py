#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from ete3 import Tree
import numpy as np
import scipy as scipy
import scipy.sparse


microbialdata=np.loadtxt("/Users/Evan/Desktop/pythontrees/microbialmatall.csv",delimiter = ",")
haarlike=scipy.sparse.load_npz("/Users/Evan/Desktop/pythontrees/97haarlike.npz")
t = Tree("/Users/Evan/Desktop/pythontrees/97_otus_unannotated.tree",format=1)
allleaves=t.get_leaves()
leafcount=len(t)
abundvec=np.zeros((leafcount,microbialdata.shape[1]-1))
node2leaves = t.get_cached_content()
N=microbialdata.shape[1]-1



#Index leaves of the tree to find OTUs
leaflist=[]
for i in range(leafcount):
    leaf=allleaves[i]
    leaflist.append(leaf.name)

#match OTUs to leaves in the tree    
for i in range(microbialdata.shape[0]):
    store=microbialdata[i,0]
    num=store.astype(int)
    name=num.astype(str)
    for j in range(1,N+1):
        abundvec[leaflist.index(name),j-1]=microbialdata[i,j]

#normalize abundance counts        
abundvecnormalized=abundvec/np.transpose(np.sum(abundvec, axis=0)[:,None])
abunds=scipy.sparse.csr_matrix(abundvecnormalized)

#Project onto the Haar-like basis and save results
mags=np.transpose(scipy.sparse.csr_matrix.transpose(abunds)@scipy.sparse.csr_matrix.transpose(haarlike))
scipy.sparse.save_npz("LOCALDIRECTORY", mags, compressed=True) 

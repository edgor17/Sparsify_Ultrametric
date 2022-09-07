#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 14:31:46 2022

@author: Evan
"""
import numpy as np
import scipy as scipy
import scipy.sparse
from skbio.stats.ordination import pcoa

mags=scipy.sparse.load_npz("/Users/Evan/Desktop/pythontrees/mags.npz")
mags=mags.todense()
N=np.shape(mags)[1]

pseudodiag=scipy.sparse.load_npz("/Users/Evan/Desktop/pythontrees/97pseudodiag.npz")
eigest=pseudodiag.diagonal()

#Build Haar-like distance matrix
D=np.zeros((N,N))
for i in range(N):
    for j in range(i+1,N):
        distdiff=((mags[:,i]-mags[:,j]))
        d=np.sum(np.multiply(np.transpose(np.square(distdiff)),eigest))
        D[i,j]=np.sqrt(d)        
D=D+np.transpose(D)
pcoa(D, method='eigh', number_of_dimensions=0, inplace=False)
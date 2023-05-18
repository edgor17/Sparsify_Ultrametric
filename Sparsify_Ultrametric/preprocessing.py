import pandas as pd
import scipy as scipy
import scipy.sparse
import numpy as np
from ete3 import Tree


def PreProcess(featuretable,metadata,labelname,labeltype,tree,haarlike):
    '''
    

    Parameters
    ----------
    featuretable : DataFrame
        OTU table
    metadata : DataFrame
        Metadata
    labelname : str
        Catagory of metadata to use
    labeltype : str
        classification or regression
    tree : TreeNode 
        Tree
    haarlike : scipy.sparse.csr_matrix
        Haar-like basis 

    Returns
    -------
    X : pd.DataFrame
        OTU abundance table
    Y : pd.Series
        Sample labels
    mags : scipy.sparse.csr_matrix
        Haar-like mags

    '''
    
    def Match_to_tree(data, tree):
        '''
        Parameters
        ----------
        data : np.ndarray
            OTU abundanes. Must include OTU ids as columns
        tree : str
            Filepath to a .nwk file

        Returns
        -------
        abundvec : np.ndarray
            OTU abundances as mapped onto the tree

        '''
        """
        Maps the given OTU abundances onto a reference tree. Data must include OTU ids as columns 
        Tree should be a filepath to a .nwk file
        """
        
        allleaves=tree.get_leaves()
        leafcount=len(tree)
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
        return abundvec

    
    featuretable=featuretable.transpose()
    X=featuretable.iloc[1: , :]
    samplenames=X.axes[0].tolist()
    
    if labeltype=='classification':
        labels=list(metadata[labelname].unique())
        dic = dict(zip(labels, list(range(1,len(labels)+1))))
        tosort=[]
        for name in samplenames:
            try:
                tosort.append(dic[metadata.loc[metadata['sample_name'] == name][labelname].values[0]])  
            except KeyError:
                try:
                    tosort.append(dic[metadata.loc[metadata['#SampleID'] == name][labelname].values[0]])  
                except:
                    tosort.append('delete')  
        X['labels']=tosort
        X=X[X.labels != 'delete']
        X=X.sort_values('labels')       
        Y=X['labels']
        X.drop(X.columns[len(X.columns)-1], axis=1, inplace=True)         
                
    elif labeltype=='regression':
        dic=None
        tosort=[]
        for name in samplenames:
            try:
                tosort.append(metadata.loc[metadata['sample_name'] == name][labelname].tolist()[0])
            except KeyError:
                tosort.append(metadata.loc[metadata['#SampleID'] == name][labelname].tolist()[0])
        floattosort=[]    
        for item in tosort:
            try:
                floattosort.append(float(item))
            except:
                floattosort.append('delete')
        X['labels']=floattosort
        X=X[X.labels != 'delete']
        X=X.sort_values('labels')       
        Y=X['labels']
        X.drop(X.columns[len(X.columns)-1], axis=1, inplace=True)    
    
    print('Mapping to Tree')
    values=X.values   
    names=list(featuretable.iloc[0,:].values)
    names=np.array(names)
    names=names.astype(int)
    look=np.vstack((names,values))
    abundvecs=Match_to_tree(look,tree)
    abunds=scipy.sparse.csr_matrix(abundvecs.T)
    
    print('Projecting onto Haar Basis')
    abunds=abunds/abunds.sum(0)
    abunds=scipy.sparse.csr_matrix(abunds)
    mags=haarlike@abunds
    
    return X,Y,mags,dic

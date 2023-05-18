import pandas as pd
import scipy as scipy
import scipy.sparse
import numpy as np
from ete3 import Tree

def find_matching_index(list1, list2):

    inverse_index = { element: index for index, element in enumerate(list1) }

    return [(index)
        for index, element in enumerate(list2) if element in inverse_index]
       

def Haar_Build(tree):
    '''
    

    Parameters
    ----------
    tree : str
        filepath to .nwk tree, in current version must be bifurcating

    Returns
    -------
    sparsehaarlike : scipy.sparse.csr_matrix
        Haar-like basis

    '''
    t = tree
    node2leaves = t.get_cached_content() #return dictionary of node instances, allows quick access to node attributes without traversing the tree
    numleaves=len(t) 
    sparsehaarlike=scipy.sparse.lil_matrix((numleaves,numleaves)) #Initialize row-based list of lists sparse matrix to store Haar-like vectors
    allleaves=t.get_leaves()
    mastersplit=t.children
    lilmat=scipy.sparse.lil_matrix((numleaves,numleaves)) #this is where we collect lstar vectors
    
    
    i=0 #ordering of nodes in post order traversal
    for node in t.traverse("postorder"):
        node.add_features(pos=find_matching_index(node,allleaves)) #store indices of leaves under each internal node
        veclen=len(node2leaves[node])
        if not node.is_leaf():
            node.add_features(loc=i) #add node index to node features
            if veclen==2:
                child=node.children
                lstar=np.zeros((numleaves,1))
                index0=child[0].pos
                index1=child[1].pos
                lstar[index0]=1
                lstar[index1]=1
                lilmat[i]=np.transpose(lstar)
                haarvec=np.zeros((numleaves,1))
                haarvec[index0]=1/np.sqrt(2)
                haarvec[index1]=-1/np.sqrt(2)
                sparsehaarlike[i]=np.transpose(haarvec)
                print(i)
                i=i+1
            else:
                child=node.children
                if len(node2leaves[child[0]])==1:
                    lstar0=np.zeros((numleaves,1))
                    index0=child[0].pos
                    lstar0[index0]=1
                    index=child[1].loc
                    lstar1=np.transpose(lilmat[index].todense())
                    index1=child[1].pos                
                    lstar1[index1]=lstar1[index1]+len(child[1])*1
                    lilmat[i]=np.transpose(lstar0)+np.transpose(lstar1)
                    L1=np.count_nonzero(lstar1)
                    haarvec=np.zeros((numleaves,1))
                    haarvec[index0]=np.sqrt(L1/(L1+1))
                    haarvec[index1]=-np.sqrt(1/(L1*(L1+1)))
                    sparsehaarlike[i]=np.transpose(haarvec)
                    print(i)
                    i=i+1
                elif len(node2leaves[child[1]])==1:
                    lstar1=np.zeros((numleaves,1))
                    index1=child[1].pos
                    lstar1[index1]=1
                    index=child[0].loc
                    lstar0=np.transpose(lilmat[index].todense())
                    index0=child[0].pos
                    lstar0[index0]=lstar0[index0]+len(child[0])*1
                    lilmat[i]=np.transpose(lstar1)+np.transpose(lstar0)
                    L0=np.count_nonzero(lstar0)
                    haarvec=np.zeros((numleaves,1))
                    haarvec[index0]=np.sqrt(1/(L0*(L0+1)))
                    haarvec[index1]=-np.sqrt(L0/((L0+1)))
                    sparsehaarlike[i]=np.transpose(haarvec)
                    print(i)
                    i=i+1
                else:
                    index0=child[0].loc
                    index1=child[1].loc
                    lstar0=np.transpose(lilmat[index0].todense())
                    lstar1=np.transpose(lilmat[index1].todense())
                    index00=child[0].pos
                    lstar0[index00]=lstar0[index00]+len(child[0])*1
                    index11=child[1].pos
                    lstar1[index11]=lstar1[index11]+len(child[1])*1
                    lilmat[i]=np.transpose(lstar0)+np.transpose(lstar1)
                    L0=np.count_nonzero(lstar0)
                    L1=np.count_nonzero(lstar1)
                    haarvec=np.zeros((numleaves,1))
                    haarvec[index00]=np.sqrt(L1/(L0*(L0+L1)))
                    haarvec[index11]=-np.sqrt(L0/(L1*(L0+L1)))
                    sparsehaarlike[i]=np.transpose(haarvec)
                    print(i)
                    i=i+1
    
    sparsehaarlike[len(allleaves)-1]=np.repeat(1/np.sqrt(len(allleaves)),len(allleaves))
    lilmat=lilmat.tocsr()
    sparsehaarlike=sparsehaarlike.tocsr()
    
    
    return sparsehaarlike


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

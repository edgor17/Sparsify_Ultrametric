import numpy as np
import scipy as scipy

#find_matching_index() was written by user Olivier Melançon in https://stackoverflow.com/questions/49247506/how-to-efficiently-find-the-indices-of-matching-elements-in-two-lists

def find_matching_index(list1, list2):

    inverse_index = { element: index for index, element in enumerate(list1) }

    return [(index)
        for index, element in enumerate(list2) if element in inverse_index]
       
def sparsify(tree):
  '''

    Parameters
    ----------
    tree : TreeNode
        Tree

    Returns
    -------
    sparsehaarlike : scipy.sparse.csr_matrix
        Haar-like wavelet basis vectors
    pseudodiag : scipy.sparse.coo_matrix
        pseudodiagonalized ultrametric matrix corresponding to the tree

    '''

  t = tree 
  node2leaves = t.get_cached_content() #return dictionary of node instances, allows quick access to node attributes without traversing the tree
  numleaves=len(t) 
  sparsehaarlike=scipy.sparse.lil_matrix((numleaves,numleaves)) #Initialize row-based list of lists sparse matrix to store Haar-like vectors
  allleaves=t.get_leaves()
  mastersplit=t.children
  lilmat=scipy.sparse.lil_matrix((numleaves,numleaves))


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
              lstar[index0]=child[0].dist
              lstar[index1]=child[1].dist
              lilmat[i]=np.transpose(lstar)
              haarvec=np.zeros((numleaves,1))
              haarvec[index0]=1/np.sqrt(2)
              haarvec[index1]=-1/np.sqrt(2)
              sparsehaarlike[i]=np.transpose(haarvec)
              i=i+1
          else:
              child=node.children
              if len(node2leaves[child[0]])==1:
                  lstar0=np.zeros((numleaves,1))
                  index0=child[0].pos
                  lstar0[index0]=child[0].dist
                  index=child[1].loc
                  lstar1=np.transpose(lilmat[index].todense())
                  index1=child[1].pos                
                  lstar1[index1]=lstar1[index1]+len(child[1])*child[1].dist
                  lilmat[i]=np.transpose(lstar0)+np.transpose(lstar1)
                  L1=np.count_nonzero(lstar1)
                  haarvec=np.zeros((numleaves,1))
                  haarvec[index0]=np.sqrt(L1/(L1+1))
                  haarvec[index1]=-np.sqrt(1/(L1*(L1+1)))
                  sparsehaarlike[i]=np.transpose(haarvec)
                  i=i+1
              elif len(node2leaves[child[1]])==1:
                  lstar1=np.zeros((numleaves,1))
                  index1=child[1].pos
                  lstar1[index1]=child[1].dist
                  index=child[0].loc
                  lstar0=np.transpose(lilmat[index].todense())
                  index0=child[0].pos
                  lstar0[index0]=lstar0[index0]+len(child[0])*child[0].dist
                  lilmat[i]=np.transpose(lstar1)+np.transpose(lstar0)
                  L0=np.count_nonzero(lstar0)
                  haarvec=np.zeros((numleaves,1))
                  haarvec[index0]=np.sqrt(1/(L0*(L0+1)))
                  haarvec[index1]=-np.sqrt(L0/((L0+1)))
                  sparsehaarlike[i]=np.transpose(haarvec)
                  i=i+1
              else:
                  index0=child[0].loc
                  index1=child[1].loc
                  lstar0=np.transpose(lilmat[index0].todense())
                  lstar1=np.transpose(lilmat[index1].todense())
                  index00=child[0].pos
                  lstar0[index00]=lstar0[index00]+len(child[0])*child[0].dist
                  index11=child[1].pos
                  lstar1[index11]=lstar1[index11]+len(child[1])*child[1].dist
                  lilmat[i]=np.transpose(lstar0)+np.transpose(lstar1)
                  L0=np.count_nonzero(lstar0)
                  L1=np.count_nonzero(lstar1)
                  haarvec=np.zeros((numleaves,1))
                  haarvec[index00]=np.sqrt(L1/(L0*(L0+L1)))
                  haarvec[index11]=-np.sqrt(L0/(L1*(L0+L1)))
                  sparsehaarlike[i]=np.transpose(haarvec)
                  i=i+1

  #Because phylogenetic trees are two ORB-trees, we handle each out rooting seperately                
  mastervec0=np.hstack((np.repeat(1/np.sqrt(len(mastersplit[0])),len(mastersplit[0])), np.repeat(0,len(mastersplit[1]))))         
  mastervec1=np.hstack((np.repeat(0,len(mastersplit[0])),np.repeat(1/np.sqrt(len(mastersplit[1])),len(mastersplit[1]))))
  sparsehaarlike[-2]=mastervec0
  sparsehaarlike[-1]=mastervec1
  side0=np.transpose(lilmat[-2,0:len(child[0])].todense())
  side1=np.transpose(lilmat[-2,len(child[0]):numleaves].todense())
  lilmat[-1]=scipy.copy(lilmat[-2].todense())
  lilmat[-2,len(child[0]):-1]=0
  lilmat[-2,-1]=0
  lilmat[-1,0:len(child[0])]=0

  lilmat=lilmat.tocsr()
  sparsehaarlike=sparsehaarlike.tocsr()

  row1=np.zeros(len(allleaves)*len(allleaves))    
  col=np.zeros(len(allleaves)*len(allleaves))
  data=np.zeros(len(allleaves)*len(allleaves))

  #scipy.sparse.save_npz("/Users/Evan/Desktop/pythontrees/97lstar", lilmat, compressed=True)


  i=0
  for node in t.traverse("postorder"):
      if not node.is_leaf():
          if not node==t:
              tempnode=node
              index=tempnode.loc
              print(index)
              row1[i]=index
              col[i]=index
              compvec=np.transpose(lilmat[index].todense())
              mulmat=sparsehaarlike[index].todense()
              #data[i]=scipy.sparse.csr_matrix.dot(scipy.sparse.csr_matrix.multiply(sparsepseudohaar[index],sparsepseudohaar[index]),np.transpose(lilmat[index])).todense()
              data[i]=np.dot(np.multiply(mulmat,mulmat),compvec)
              i=i+1
              while not tempnode.up==t:
                  tempnode=tempnode.up
                  index2=tempnode.loc
                  row1[i]=index2
                  col[i]=index
                  #data[i]=scipy.sparse.csr_matrix.dot(scipy.sparse.csr_matrix.multiply(sparsepseudohaar[index2],sparsepseudohaar[index]),np.transpose(lilmat[index])).todense()
                  data[i]=np.dot(np.multiply(sparsehaarlike[index2].todense(),mulmat),compvec)
                  i=i+1
  i=i+1
  for node in t.traverse("postorder"):
      if not node.is_leaf():    
          if not node==t:         
              index2=node.loc
              row1[i]=numleaves-2
              col[i]=index2          
              data[i]=np.dot(np.multiply(sparsehaarlike[index2].todense(),sparsehaarlike[-2].todense()),np.transpose(lilmat[-2].todense()))
              i=i+1
          if node==t:
              row1[i]=numleaves-2
              col[i]=numleaves-2
              data[i]=np.dot(np.multiply(sparsehaarlike[-2].todense(),sparsehaarlike[-2].todense()),np.transpose(lilmat[-2].todense()))
              i=i+1
  i=i+1
  for node in t.traverse("postorder"):
      if not node.is_leaf():    
          if not node==t:         
              index2=node.loc
              row1[i]=numleaves-1
              col[i]=index2
              data[i]=np.dot(np.multiply(sparsehaarlike[index2].todense(),sparsehaarlike[-1].todense()),np.transpose(lilmat[-1].todense()))
              i=i+1
          if node==t:
              row1[i]=numleaves-1
              col[i]=numleaves-1
              data[i]=np.dot(np.multiply(sparsehaarlike[-1].todense(),sparsehaarlike[-1].todense()),np.transpose(lilmat[-1].todense()))
              i=i+1

  row1=row1[0:i]
  col=col[0:i]
  data=data[0:i]       

  a=scipy.sparse.coo_matrix((data, (row1, col)), shape=(numleaves, numleaves))




  b=a.diagonal()

  rowdiag=np.zeros(len(allleaves))    
  coldiag=np.zeros(len(allleaves))
  datadiag=np.zeros(len(allleaves))

  for j in range(len(allleaves)):
      rowdiag[j]=j
      coldiag[j]=j
      datadiag[j]=b[j]

  c=scipy.sparse.coo_matrix((datadiag, (rowdiag, coldiag)), shape=(len(allleaves), len(allleaves)))        
  pseudodiag=a+a.transpose()-c
  
  return sparsehaarlike, pseudodiag

#Run this to load data prior to producing any plots

from Sparsify_Ultrametric import Sparsify
from Sparsify_Ultrametric.utils import compute_Haar_dist, PCoA
from Sparsify_Ultrametric.preprocessing import PreProcess
from ete3 import Tree
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.pyplot import cm
import numpy as np

tree = Tree("Sparsify_Ultrametric/raw_data/97_otus_unannotated.tree",format=1)
haarlike=scipy.sparse.load_npz('Sparsify_Ultrametric/precomputed/97haarlike.npz')
pseudodiag=scipy.sparse.load_npz('Sparsify_Ultrametric/precomputed/97pseudodiag.npz')
featuretable=pd.read_csv("Sparsify_Ultrametric/raw_data/otus.txt", sep='\t')
metadata=pd.read_csv("Sparsify_Ultrametric/raw_data/metadata.txt", sep='\t')
[X,Y,mags,dic]=PreProcess(featuretable,metadata,'end_depth','regression',tree,haarlike)
lambdav=scipy.sparse.csr_matrix.diagonal(pseudodiag)
[D,modmags]=compute_Haar_dist(mags,lambdav)






#Figure 5(a)

m = scipy.sparse.coo_matrix(pseudodiag)
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
ax.plot(m.col, np.zeros(len(m.col)), 's', color='#1f77b4', marker=',',linewidth=0)
ax.set_xlim(0, m.shape[1])
ax.set_ylim(0, m.shape[0])
ax.invert_yaxis()
ax.add_patch(
     patches.Rectangle(
        (0, 0),
        2446,
        2446,
        fill=True      
     ) ) 
ax.add_patch(
     patches.Rectangle(
        (2446, 2446),
        96876,
        96876,
        fill=True      
     ) ) 
for spine in ax.spines.values():
  spine.set_visible(False)
ax.tick_params(axis=u'both', which=u'both',length=0)


#Figure 5(b)

fig = plt.figure()
ax = fig.add_subplot(111, facecolor='white')
ax.plot(m.col, m.row, 's', color='#1f77b4', marker=',',linewidth=0)
ax.set_xlim(0, m.shape[1])
ax.set_ylim(0, m.shape[0])
ax.set_aspect('equal')
for spine in ax.spines.values():
  spine.set_visible(False)
ax.invert_yaxis()
ax.set_aspect('equal')
#ax.set_xticks([])
#ax.set_yticks([])
ax.tick_params(axis=u'both', which=u'both',length=0)
ax.figure.show()

#Produces figure 6

eigvals=scipy.sparse.linalg.eigsh(pseudodiag,500,return_eigenvectors=False)
fig,ax = plt.subplots()
true=ax.scatter(np.linspace(1,500,500),np.log10(np.flip(eigvals)),facecolors='none',s=20,edgecolors='#1f77b4',label='True Eigenvalues')
approx=ax.scatter(np.linspace(1,500,500),np.log10(np.flip(np.sort(lambdav))[0:500]),s=20,marker='x',c='#ff7f0e',label='Similar Matrix Diagonal')
plt.rcParams["mathtext.fontset"] = "cm"
plt.xlabel('Eigenvalue Rank')
plt.ylabel(r'$\log_{10} \lambda$')
ax.legend(handles=[true,approx],loc=1,prop={'size': 7.5})
axins = ax.inset_axes([0.15, 0.4, 0.47, 0.47])
axins.scatter(np.linspace(1,500,500),np.log10(np.flip(eigvals)),facecolors='none',s=10,edgecolors='#1f77b4',label='True Eigenvalues')
axins.scatter(np.linspace(1,500,500),np.log10(np.flip(np.sort(lambdav))[0:500]),s=10,marker='x',c='#ff7f0e',label='Similar Matrix Diagonal')
x1, x2, y1, y2 = -0, 30, 2, 5.2
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.xaxis.set_ticks_position('none') 
axins.yaxis.set_ticks_position('none') 
axins.set_xticklabels([])
axins.set_yticklabels([])
ax.indicate_inset_zoom(axins, edgecolor="grey",linestyle='dashed')
axins.spines['top'].set_edgecolor('grey')
axins.spines['right'].set_edgecolor('grey')
axins.spines['left'].set_edgecolor('grey')
axins.spines['bottom'].set_edgecolor('grey')
axins.spines['top'].set_linestyle('dashed')
axins.spines['right'].set_linestyle('dashed')
axins.spines['left'].set_linestyle('dashed')
axins.spines['bottom'].set_linestyle('dashed')
plt.title('500 Largest Eigenvalues of 97% Greengenes')



#Figure 7

unifrac=pd.read_csv("Sparsify_Ultrametric/precomputed/unifracdists.csv", sep=',')
DPCoA=pd.read_csv("Sparsify_Ultrametric/precomputed/DPCoAdists.csv", sep=',')
depths=[]
for sample in unifrac.columns:
    depths.append(metadata.loc[metadata['sample_name'] == sample]['end_depth'].tolist()[0])
depths=np.array(depths)
indices=np.argsort(depths)
depths=depths[indices]
sortedunifrac=unifrac.values[np.ix_(indices,indices)]
sortedDPCoA=DPCoA.values[np.ix_(indices,indices)]
unifraccoord=PCoA(sortedunifrac,2)
DPCoAcoord=PCoA(sortedDPCoA,2)
haarlikecoord=PCoA(D,2)
color = iter(cm.viridis(np.linspace(0, 1, 9)))
fig,ax=plt.subplots(3,figsize=(8, 8))
for item in np.unique(depths):
    c=next(color)
    unifracplot=ax[0].scatter(unifraccoord[np.where(depths==item),0] ,unifraccoord[np.where(depths==item),1],c=c,label='.00'+str(item))
    DPCoAplot=ax[1].scatter(-DPCoAcoord[np.where(depths==item),0] ,DPCoAcoord[np.where(depths==item),1],c=c,label='.00'+str(item))
    haarlikeplot=ax[2].scatter(-haarlikecoord[np.where(depths==item),0] ,haarlikecoord[np.where(depths==item),1],c=c,label='.00'+str(item))
ax[0].title.set_text('UniFrac Embedding')
ax[1].title.set_text('DPCoA Embedding')
ax[2].title.set_text('Haar-like Distance Embedding')
plt.legend(title="Depth in Meters",loc='center left', bbox_to_anchor=(1.1, 1.85))
plt.subplots_adjust(hspace=.4)


#Produces figure 8(a)

diff=np.square(modmags[0,:].todense()-modmags[16,:].todense()).T
diff[diff==0]=np.nan
fig, ax = plt.subplots(figsize=(10,5))
markerline, stemline, baseline, = ax.stem(np.linspace(1,99322,99322), diff,linefmt='k-',markerfmt='ko',basefmt='k.')
plt.setp(stemline, linewidth = 1.25)
plt.setp(stemline, 'linestyle', 'dotted')
plt.setp(markerline, markersize = 2, fillstyle='none')
one=plt.scatter(99310,diff[99310].item(),s=50,marker='D',c='#2a4858')
two=plt.scatter(6078,diff[6078].item(),s=50,marker='D',c='#00898a')
three=plt.scatter(67316,diff[67316].item(),s=50,marker='D',c='#64c987')
plt.rcParams["mathtext.fontset"] = "cm"
plt.xlabel('Internal Node Index (Postorder Traversal)')
plt.ylabel(r'$\lambda_v \Delta_v^2$')
plt.title('Comparing the Shallowest and Deepest Mat Samples via the Haar-like Basis')
plt.legend([one,two,three],[r'$\varphi_{99311}\sim 4.84\times10^{-2}$',r'$\varphi_{6079}\sim 4.84\times10^{-3}$',r'$\varphi_{67317}\sim 4.65\times10^{-3}$'])
plt.show()

#Figure 8(b) was produced in iTOL. To reproduce this figure, upload the .nwk file provided in /Precomputed. Next drag all of the .txt files in /Precomputed onto the tree visualization

diff=np.square(modmags[0,:].todense()-modmags[16,:].todense()).T
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
plt.legend([one,two,three],[r'$\varphi_{99311}\sim 4.84\times10^{-2}$',r'$\varphi_{6079}\sim 4.84\times10^{-3}$',r'$\varphi_{67317}\sim 4.75\times10^{-3}$'])
plt.savefig('/Users/Evan/Desktop/pythontrees/paper1spectrogram.png',dpi=400)

plt.show()


fig,ax = plt.subplots()
true=ax.scatter(np.linspace(1,500,500),np.log10((eigvals[0])),facecolors='none',s=20,edgecolors='#1f77b4',label='True Eigenvalues')
approx=ax.scatter(np.linspace(1,500,500),np.log10(np.flip(np.sort(lambdav))[0:500]),s=20,marker='x',c='#ff7f0e',label='Similar Matrix Diagonal')
plt.rcParams["mathtext.fontset"] = "cm"
plt.xlabel('Eigenvalue Rank')
plt.ylabel(r'$\log_{10} \lambda$')
ax.legend(handles=[true,approx],loc=1,prop={'size': 7.5})
axins = ax.inset_axes([0.15, 0.4, 0.47, 0.47])
axins.scatter(np.linspace(1,500,500),np.log10((eigvals[0])),facecolors='none',s=10,edgecolors='#1f77b4',label='True Eigenvalues')
axins.scatter(np.linspace(1,500,500),np.log10(np.flip(np.sort(lambdav))[0:500]),s=10,marker='x',c='#ff7f0e',label='Similar Matrix Diagonal')
x1, x2, y1, y2 = -0, 30, 2, 5.2
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.xaxis.set_ticks_position('none') 
axins.yaxis.set_ticks_position('none') 
axins.set_xticklabels([])
axins.set_yticklabels([])
ax.indicate_inset_zoom(axins, edgecolor="black")
plt.title('500 Largest Eigenvalues of 97% Greengenes')
plt.savefig('/Users/Evan/Desktop/pythontrees/paper1eigapprox.png',dpi=400)


from matplotlib.pyplot import cm


metadata=pd.read_csv("/Users/Evan/Desktop/mat/metadata.txt", sep='\t')
unifrac=pd.read_csv("/Users/Evan/Desktop/pythontrees/unifracdists.csv", sep=',')
DPCoA=pd.read_csv("/Users/Evan/Desktop/pythontrees/DPCoAdists.csv", sep=',')
depths=[]
for sample in unifrac.columns:
    depths.append(metadata.loc[metadata['sample_name'] == sample]['end_depth'].tolist()[0])
depths=np.array(depths)
indices=np.argsort(depths)
depths=depths[indices]
sortedunifrac=unifrac.values[np.ix_(indices,indices)]
sortedDPCoA=DPCoA.values[np.ix_(indices,indices)]

featuretable=pd.read_csv("/Users/Evan/Desktop/mat/otus.txt", sep='\t')
metadata=pd.read_csv("/Users/Evan/Desktop/mat/metadata.txt", sep='\t')
[X,Y,mags,dic]=PreProcess(featuretable,metadata,'end_depth','regression',tree,haarlike)
lambdav=scipy.sparse.csr_matrix.diagonal(pseudodiag)
[haarlikedist,modmags]=compute_Haar_dist(mags,lambdav)

unifraccoord=PCoA(sortedunifrac,2)
DPCoAcoord=PCoA(sortedDPCoA,2)
haarlikecoord=PCoA(haarlikedist,2)

color = iter(cm.viridis(np.linspace(0, 1, 9)))
fig,ax=plt.subplots(3,figsize=(8, 8))
for item in np.unique(depths):
    c=next(color)
    unifracplot=ax[0].scatter(unifraccoord[np.where(depths==item),0] ,-unifraccoord[np.where(depths==item),1],c=c,label='.00'+str(item))
    DPCoAplot=ax[1].scatter(-DPCoAcoord[np.where(depths==item),0] ,-DPCoAcoord[np.where(depths==item),1],c=c,label='.00'+str(item))
    haarlikeplot=ax[2].scatter(-haarlikecoord[np.where(depths==item),0] ,-haarlikecoord[np.where(depths==item),1],c=c,label='.00'+str(item))
ax[0].title.set_text('UniFrac Embedding')
ax[1].title.set_text('DPCoA Embedding')
ax[2].title.set_text('Haar-like Distance Embedding')
plt.legend(title="Depth in Meters",loc='center left', bbox_to_anchor=(1.1, 1.85))
plt.subplots_adjust(hspace=.4)
plt.savefig('/Users/Evan/Desktop/pythontrees/paper1embeddings.png',dpi=400, bbox_inches='tight')





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
plt.savefig('/Users/Evan/Desktop/pythontrees/paper1dense.png',dpi=400)




pseudodiag=scipy.sparse.load_npz('/precomputed/97pseudodiag.npz')
m = coo_matrix(pseudodiag)
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
plt.savefig('/Users/Evan/Desktop/pythontrees/paper1sparsity.png',dpi=400)

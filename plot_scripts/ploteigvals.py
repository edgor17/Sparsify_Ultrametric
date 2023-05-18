#To generate figure 7

pseudodiag=scipy.sparse.load_npz('/precomputed/97pseudodiag.npz')
eigvals=scipy.sparse.linalg.eigs(pseudodiag,500)

fig,ax = plt.subplots()
true=ax.scatter(np.linspace(1,500,500),np.log(np.flip(eigvals)),label='True Eigenvalues')
approx=ax.scatter(np.linspace(1,500,500),np.log(np.flip(np.sort(lambdav))[0:500]),label='Similar Matrix Diagonal')
plt.xlabel('Eigenvalue Rank')
ax.legend(handles=[true,approx])
plt.title('500 Largest Eigenvalues of 97% Greengenes')

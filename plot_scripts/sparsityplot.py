#used to generate figure 6

pseudodiag=scipy.sparse.load_npz('97pseudodiag.npz')
m = coo_matrix(pseudodiag)
fig = plt.figure()
ax = fig.add_subplot(111, facecolor='white')
ax.plot(m.col, m.row, 's', color='black', marker=',',linewidth=0)
ax.set_xlim(0, m.shape[1])
ax.set_ylim(0, m.shape[0])
ax.set_aspect('equal')
for spine in ax.spines.values():
  spine.set_visible(False)
ax.invert_yaxis()
ax.set_aspect('equal')
ax.set_xticks([])
ax.set_yticks([])
ax.figure.show()

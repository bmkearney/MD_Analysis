import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import MDAnalysis as mda
#from MDAnalysis.tests.datafiles import PSF, DCD
from MDAnalysis.analysis import rms, pca, align

import warnings
# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

#u = mda.Universe(PSF, DCD)
PDB="step3_input.pdb"
DCD="step5_1.dcd"

u=mda.Universe(PDB,DCD)

#We can generate an average structure to align to with the align.AverageStructure class.
#Here we first align to the first frame (ref_frame=0), and then average the coordinates.
average = align.AverageStructure(u, u, select='protein and name CA',
                                 ref_frame=0).run()
ref = average.results.universe

aligner = align.AlignTraj(u, ref,
                          select='protein and name CA',
                          in_memory=True).run()

c_alphas = u.select_atoms('protein and name CA')
R = rms.RMSF(c_alphas).run()
protein = u.select_atoms("protein")

plt.plot(c_alphas.resids, R.results.rmsf)
plt.xlabel('Residue number NL')
plt.ylabel('RMSF ($\AA$)')
#plt.axvspan(122, 159, zorder=0, alpha=0.2, color='orange', label='LID')
#plt.axvspan(30, 59, zorder=0, alpha=0.2, color='green', label='NMP')
plt.axvspan(29, 50, zorder=0, alpha=0.2, color='green', label='Helix 1')
plt.axvspan(53, 76, zorder=0, alpha=0.2, color='orange', label='Helix 2')
plt.axvspan(81, 105, zorder=0, alpha=0.2, color='green', label='Helix 3')
plt.axvspan(113, 130, zorder=0, alpha=0.2, color='orange', label='Helix 4')
plt.axvspan(137, 150, zorder=0, alpha=0.2, color='green', label='Helix 5')
plt.axvspan(153, 160, zorder=0, alpha=0.2, color='orange', label='Helix 6')
plt.axvspan(165, 180, zorder=0, alpha=0.2, color='green', label='Helix 7')
ax = plt.gca()
ax.set_ylim([0, 4])
plt.legend()
plt.savefig('RMSF.png',dpi=1200)
plt.show()

with mda.Writer("align.dcd", protein.n_atoms) as W:
    for ts in u.trajectory:
        W.write(protein)

with mda.Writer("align.pdb", protein.n_atoms) as W:
        W.write(protein)



"""
aligner = align.AlignTraj(u, u, select='backbone',
                          in_memory=True).run()

pc = pca.PCA(u, select='backbone',
             align=True, mean=None,
             n_components=None).run()

backbone = u.select_atoms('backbone')
n_bb = len(backbone)
print('There are {} backbone atoms in the analysis'.format(n_bb))
print(pc.p_components.shape)


pc.variance[0]

for i in range(3):
    print(f"Cumulated variance: {pc.cumulated_variance[i]:.3f}")
plt.plot(pc.cumulated_variance[:10])
plt.xlabel('Principal component')
plt.ylabel('Cumulative variance')
plt.show()

transformed = pc.transform(backbone, n_components=3)
transformed.shape

df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(3)])
df['Time (ps)'] = df.index * u.trajectory.dt
df.head()
"""



import MDAnalysis as mda
from MDAnalysis.analysis.base import (AnalysisBase,
                                      AnalysisFromFunction,
                                      analysis_class)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def radgyr(atomgroup, masses, total_mass=None):
    # coordinates change for each frame
    coordinates = atomgroup.positions
    center_of_mass = atomgroup.center_of_mass()

    # get squared distance from center
    ri_sq = (coordinates-center_of_mass)**2
    # sum the unweighted positions
    sq = np.sum(ri_sq, axis=1)
    sq_x = np.sum(ri_sq[:,[1,2]], axis=1) # sum over y and z
    sq_y = np.sum(ri_sq[:,[0,2]], axis=1) # sum over x and z
    sq_z = np.sum(ri_sq[:,[0,1]], axis=1) # sum over x and y

    # make into array
    sq_rs = np.array([sq, sq_x, sq_y, sq_z])

    # weight positions
    rog_sq = np.sum(masses*sq_rs, axis=1)/total_mass
    # square root and return
    return np.sqrt(rog_sq)

PDB="step3_input.pdb"
DCD="step5_1.dcd"

u=mda.Universe(PDB,DCD)
protein = u.select_atoms('protein')

rog = AnalysisFromFunction(radgyr, u.trajectory,
                           protein, protein.masses,
                           total_mass=np.sum(protein.masses))
rog.run();
rog.results['timeseries'].shape
labels = ['all', 'x-axis', 'y-axis', 'z-axis']
for col, label in zip(rog.results['timeseries'].T, labels):
    plt.plot(col, label=label)
plt.legend()
plt.ylabel('Radius of gyration (Å)')
plt.xlabel('Frame (NO Ligand)')
plt.show()

rog_10 = AnalysisFromFunction(radgyr, u.trajectory,
                              protein, protein.masses,
                              total_mass=np.sum(protein.masses))

rog_10.run(start=10, stop=800, step=1)
rog_10.results['timeseries'].shape
for col, label in zip(rog_10.results['timeseries'].T, labels):
    plt.plot(col, label=label)
plt.legend()
plt.ylabel('Radius of gyration (Å)')
plt.xlabel('Frame (NO Ligand)');
plt.show()

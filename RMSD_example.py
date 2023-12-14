import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
#from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')

PDB="step3_input.pdb"
DCD="step5_1.dcd"
#u1=mda.Universe("TIPE1_no_ligand.pdb", "TIPE1_no_ligand.dcd")
mobile=mda.Universe(PDB,DCD)
ref=mda.Universe(PDB,DCD)
mobile.trajectory[-1]  # set mobile trajectory to last frame
ref.trajectory[0]  # set reference trajectory to first frame

mobile_ca = mobile.select_atoms('name CA')
ref_ca = ref.select_atoms('name CA')
unaligned_rmsd = rms.rmsd(mobile_ca.positions, ref_ca.positions, superposition=False)
print(f"Unaligned RMSD: {unaligned_rmsd:.2f}")

aligner = align.AlignTraj(mobile, ref, select='name CA', in_memory=True).run()
mobile.trajectory[-1]  # set mobile trajectory to last frame
ref.trajectory[0]  # set reference trajectory to first frame

mobile_ca = mobile.select_atoms('name CA')
ref_ca = ref.select_atoms('name CA')
aligned_rmsd = rms.rmsd(mobile_ca.positions, ref_ca.positions, superposition=False)

print(f"Aligned RMSD: {aligned_rmsd:.2f}")

with mda.Writer(pdbtrj, multiframe=True, bonds=None, n_atoms=u.atoms.n_atoms) as PDB2:
    print("Hi")

#u2=mda.Universe("TIPE1_no_ligand.dcd")
#print (u2)
#print(mda.Universe(PSF, DCD))
#print("Using MDAnalysis version", mda.__version__)

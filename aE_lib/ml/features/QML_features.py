# define functions to produce QML based features
import qml
from util.flag_handler.hdl_targetflag import flag_to_target
import numpy as np


# all features for QML are of the shape: [n, p, r]
# where n is the number of "features", each n is a single NMR environment
# p is the number of sub environment features needed in each NMR environment
# r is the representation itself

def get_FCHL_features(mols, targetflag='CCS', cutoff=5.0, max=400):

	target = flag_to_target(targetflag)

	x = []
	y = []
	r = []

	for mol in mols:
		reps = qml.fchl.generate_representation(mol.xyz, mol.types, max, cut_distance=cutoff)

		if len(target) == 1:
			for i in range(len(mol.types)):
				if mol.types[i] == target[0]:
					x.append([reps[i]])
					y.append(mol.shift[i])
					r.append([mol.molid, i])

		if len(target) == 3:
			for i in range(len(mol.types)):
				for j in range(len(mol.types)):
					if i == j:
						continue
					if mol.types[i] != target[1] and mol.types[j] != target[2]:
						continue

					if mol.coupling_len[i][j] != target[0]:
						continue

					x.append([reps[i], reps[j]])
					y.append(mol.coupling[i][j])
					r.append([mol.molid, i, j])

	return x, y, r

def get_aSLATM_mbtypes(mols):

	if len(mols) == 0:
		mbtypes = [[1],[1,1], [1,1,1], [1,1,6], [1,1,7], [1,1,8], [1,1,9], [1,6], [1,6,1], [1,6,6], [1,6,7], [1,6,8], [1,6,9], [1,7], [1,7,1], [1,7,6], [1,7,7], [1,7,8], [1,7,9], [1,8], [1,8,1], [1,8,6], [1,8,7], [1,8,8], [1,8,9], [1,9], [1,9,1], [1,9,6], [1,9,7], [1,9,8], [1,9,9], [6], [6,1], [6,1,1], [6,1,6], [6,1,7], [6,1,8], [6,1,9], [6,6], [6,6,1], [6,6,6], [6,6,7], [6,6,8], [6,6,9], [6,7], [6,7,1], [6,7,6], [6,7,7], [6,7,8], [6,7,9], [6,8], [6,8,1], [6,8,6], [6,8,7], [6,8,8], [6,8,9], [6,9], [6,9,1], [6,9,6], [6,9,7], [6,9,8], [6,9,9], [7],[7,1], [7,1,1], [7,1,6], [7,1,7], [7,1,8], [7,1,9], [7,6], [7,6,1], [7,6,6], [7,6,7], [7,6,8], [7,6,9], [7,7], [7,7,1], [7,7,6], [7,7,7], [7,7,8], [7,7,9], [7,8], [7,8,1], [7,8,6], [7,8,7], [7,8,8], [7,8,9], [7,9], [7,9,1], [7,9,6], [7,9,7], [7,9,8], [7,9,9], [8], [8,1], [8,1,1], [8,1,6], [8,1,7], [8,1,8], [8,1,9], [8,6], [8,6,1], [8,6,6], [8,6,7], [8,6,8], [8,6,9], [8,7], [8,7,1], [8,7,6], [8,7,7], [8,7,8], [8,7,9], [8,8], [8,8,1], [8,8,6], [8,8,7], [8,8,8], [8,8,9], [8,9], [8,9,1], [8,9,6], [8,9,7], [8,9,8], [8,9,9], [9], [9,1], [9,1,1], [9,1,6], [9,1,7], [9,1,8], [9,1,9], [9,6], [9,6,1], [9,6,6], [9,6,7], [9,6,8], [9,6,9], [9,7], [9,7,1], [9,7,6], [9,7,7], [9,7,8], [9,7,9], [9,8], [9,8,1], [9,8,6], [9,8,7], [9,8,8], [9,8,9], [9,9], [9,9,1], [9,9,6], [9,9,7], [9,9,8], [9,9,9]]

	else:
		nuclear_charges = []
		for tmp_mol in mols:
			nuclear_charges.append(tmp_mol.types)
		mbtypes = qml.representations.get_slatm_mbtypes(nuclear_charges)

	return mbtypes

def get_aSLATM_features(mols, targetflag='CCS', cutoff=5.0, max=400, mbtypes=[]):

		target = flag_to_target(targetflag)

		x = []
		y = []
		r = []

		if len(mbtypes) == 0:
			mbtypes = get_aSLATM_mbtypes([])


		for mol in mols:
			reps = qml.representations.generate_slatm(mol.xyz, mol.types, mbtypes,
					unit_cell=None, local=True, sigmas=[0.05,0.05], dgrids=[0.03,0.03],
					rcut=cutoff, alchemy=False, pbc='000', rpower=6)

			reps = np.asarray(reps)

			if len(target) == 1:
				for i in range(len(mol.types)):
					if mol.types[i] == target[0]:
						x.append([reps[i]])
						y.append(mol.shift[i])
						r.append([mol.molid, i])

			if len(target) == 3:
				for i in range(len(mol.types)):
					for j in range(len(mol.types)):
						if i == j:
							continue
						if mol.types[i] != target[1] and mol.types[j] != target[2]:
							continue

						if mol.coupling_len[i][j] != target[0]:
							continue

						x.append([reps[i], reps[j]])
						y.append(mol.coupling[i][j])
						r.append([mol.molid, i, j])

		return x, y, r


def get_CMAT_features(mols, targetflag='CCS', cutoff=5.0, max=400):

	target = flag_to_target(targetflag)

	x = []
	y = []
	r = []
	for mol in mols:
		if len(mol.types) > max:
			max = len(mol.types)

	for mol in mols:
		reps = qml.representations.generate_atomic_coulomb_matrix(mol.types, mol.xyz, size = max, sorting = "distance",
						central_cutoff = cutoff, central_decay = -1, interaction_cutoff = 1e6, interaction_decay = -1,
						indices = None)

		if len(target) == 1:
			for i in range(len(mol.types)):
				if mol.types[i] == int(target[0]):
					x.append([reps[i]])
					y.append(mol.shift[i])
					r.append([mol.molid, i])

		if len(target) == 3:
			for i in range(len(mol.types)):
				for j in range(len(mol.types)):
					if i == j:
						continue

					if mol.types[i] != target[1] and mol.types[j] != target[2]:
						continue

					if mol.coupling_len[i][j] != target[0]:
						continue

					x.append([reps[i], reps[j]])
					y.append(mol.coupling[i][j])
					r.append([mol.molid, i, j])

	return x, y, r

def get_ACSF_features(mols, targetflag='CCS', cutoff=5.0, max=400, elements=[]):

	target = flag_to_target(targetflag)

	x = []
	y = []
	r = []
	'''
	elements = set()
	for tmp_mol in mols:
		elements = elements.union(tmp_mol.types)
	elements = sorted(list(elements))
	'''
	for mol in mols:
		# need to put in defaults
		reps = qml.representations.generate_acsf(mol.types, mol.xyz, elements=elements)
												#, nRs2=nRs2, nRs3=nRs3,
												#nTs=nTs, eta2=eta2, eta3=eta3, zeta=zeta, rcut=cutoff, acut=acut,
												#bin_min=bin_min, gradients=gradients)

		if len(target) == 1:
			for i in range(len(mol.types)):
				if mol.types[i] == target[0]:
					x.append([reps[i]])
					y.append(mol.shift[i])
					r.append([mol.molid, i])

		if len(target) == 3:
			for i in range(len(mol.types)):
				for j in range(len(mol.types)):
					if i == j:
						continue
					if mol.types[i] != target[1] and mol.types[j] != target[2]:
						continue

					if mol.coupling_len[i][j] != target[0]:
						continue

					x.append([reps[i], reps[j]])
					y.append(mol.coupling[i][j])
					r.append([mol.molid, i, j])

	return x, y, r

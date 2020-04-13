import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../../')

from autoenrich.molecule.dataset import dataset

from . import dummy_mol as dummymol


def get_dummy_dataset(ml_size=10, at_size=10, target=[1]):

	dset = dataset()

	for i in range(ml_size):
		dset.mols.append(dummymol.get_random_mol(size=at_size))
		dset.files.append('file_', str(i))

	x = []
	y = []
	r = []

	for mol in mols:
		if len(target) == 1:
			for i in range(len(mol.types)):
				if mol.types[i] == target[0]:
					y.append(mol.shift[i])
					r.append([mol.molid, i])

		if len(target) == 3:
			for i in range(len(mol.types)):
				for j in range(len(mol.types)):
					if i == j:
						continue

					if not ( mol.types[i] == target[1] and mol.types[j] == target[2] ):
						continue

					if mol.coupling_len[i][j] != target[0]:
						continue

					y.append(mol.coupling[i][j])
					r.append([mol.molid, i, j])

	dset.x = x
	dset.y = y
	dset.r = r

	return dset



















##

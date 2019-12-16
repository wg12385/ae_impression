
from util.flag_handler.hdl_targetflag import flag_to_target
import numpy as np


def get_dummy_features(mols, targetflag='CCS'):

	target = flag_to_target(targetflag)

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

	return x, y, r

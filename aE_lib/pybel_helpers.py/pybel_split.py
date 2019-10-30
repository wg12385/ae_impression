import numpy as np
import openbabel
from reference.periodic_table import Get_periodic_table

def mol_iswhole(mol):
	iswhole = True
	atoms = mol_getatoms(mol)
	for i in range(atoms):
		for j in range(atoms):
			if i == j:
				continue
			good = False

			while good == False:
				for length in range(0, atoms):
					paths = mol_find_all_paths(mol, i, j, length)
					if len(paths) > 0:
						good = True
						break
				if not good:
					iswhole = False
					return iswhole
	return iswhole

def mol_splitmol(mol):
	atoms = mol_getatoms(mol)
	structures = []
	assigned = []
	structures.append([])

	depth = 0
	while len(assigned) < atoms:

		for i in range(atoms):
			if i in assigned:
				continue

			if len(structures[depth]) == 0:
				structures[depth].append(i)
				assigned.append(i)

			for j in range(atoms):
				if j in assigned or j == i:
					continue

				for length in range(0, atoms):
					paths = mol_find_all_paths(mol, i, j, length)
					if len(paths) > 0:
						structures[depth].append(j)
						assigned.append(j)
						break
						continue

			structures.append([])
			depth += 1

	lent = 0
	for structure in structures:
		lent += len(structure)

	assert lent == atoms

	return structures

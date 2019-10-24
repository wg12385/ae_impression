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

def mol_read_type(mol):
	type_list = []
	type_array = np.zeros(len(mol.atoms), dtype=np.int32)
	Periodic_table = Get_periodic_table()
	for i in range(len(mol.atoms)):
		type = int(mol.atoms[i].atomicnum)
		type_array[i] = type
		type_list.append(Periodic_table[type])

	return type_list, type_array

def mol_getatoms(mol):
	atomnumber = len(mol.atoms)
	return atomnumber

def mol_read_xyz(mol):
	xyz_array = np.zeros((len(mol.atoms),3), dtype=np.float64)
	for i in range(len(mol.atoms)):
		xyz_array[i][0] = float(mol.atoms[i].coords[0])
		xyz_array[i][1] = float(mol.atoms[i].coords[1])
		xyz_array[i][2] = float(mol.atoms[i].coords[2])

	# Return array is zero indexed
	return xyz_array

def mol_read_dist(mol):
	xyz_array = mol_read_xyz(mol)
	atoms = mol_getatoms(mol)
	d_array = np.zeros((atoms, atoms), dtype=np.float64)
	for i in range(atoms):
		for j in range(atoms):
			d_array[i][j] = np.absolute(np.linalg.norm(xyz_array[i] - xyz_array[j]))

	return d_array

def mol_find_all_paths(mol, start, end, coupling_length, path=[]):
		# append atom to start
		path = path + [start]
		# check if we have reached target atom
		if start == end:
			# if we have, return succesful path
			return [path]
		# define new path
		paths = []
		# loop over neighbouring atoms
		for nbr_atom in openbabel.OBAtomAtomIter(mol.atoms[start].OBAtom):
			# get ID of neighbour
			node = nbr_atom.GetId()
			# check the neighbour is not already in the path, and that the path is not over the required length
			if node not in path and len(path) <= coupling_length:
				# get new paths for the neighbour
				newpaths = mol_find_all_paths(mol, node, end, coupling_length, path)
				#for each new path, check for paths of correct length
				for newpath in newpaths:
					if len(newpath) == coupling_length+1:
						paths.append(newpath)
		return paths

def mol_get_bond_table(mol):

	atoms = len(mol.atoms)

	bond_table = np.zeros((atoms, atoms), dtype=np.int32)

	for atom1 in range(atoms):
		for atom2 in range(atom1, atoms):

			for nbr_atom in openbabel.OBAtomAtomIter(mol.atoms[atom1].OBAtom):
				check = nbr_atom.GetId()
				if atom2 != check:
					continue

				bond = mol.atoms[atom1].OBAtom.GetBond(nbr_atom)
				order = bond.GetBondOrder()

				bond_table[atom1][atom2] = int(order)
				bond_table[atom2][atom1] = int(order)

	return bond_table

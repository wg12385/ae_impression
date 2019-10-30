import numpy as np
import openbabel
from reference.periodic_table import Get_periodic_table

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

def get_coupling_lengths(mol, types, maxlen=6):
	coupling_len = np.zeros((len(types), len(types)), dtype=np.int32)
	for atom1 in range(len(types)):
		for atom2 in range(len(types)):
			coupling_paths = []
			for i in range(maxlen):
				coupling_paths.extend(mol_find_all_paths(mol, atom1, atom2, i))
			length = 999
			for path in coupling_paths:
				if len(path) < length and len(path) != 0:
					length = len(path)

			if length > maxlen+1:
				length = 0

			coupling_len[atom1][atom2] = length
			coupling_len[atom2][atom1] = length

	return coupling_len

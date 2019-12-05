# get features for TransForMer models
from reference.periodic_table import Get_periodic_table
from util.flag_handler.hdl_target_flag import flag_to_target

import ml.models.features.BCAI_calc.mol_graph_setup as BCAI

import numpy as np
import pandas as pd
import sys


def get_BCAI_features(mols, targetflag='CCS'):

	target = flag_to_target(targetflag)
	p_table = Get_periodic_table()

	# construct dataframe as BCAI requires from mols
	# atoms has: molecule_name, atom, labeled atom,
	molecule_name = [] 	# molecule name
	atom_index = []		# atom index
	atom = []			# atom type (letter)
	x = []				# x coordinate
	y = []				# y coordinate
	z = []				# z coordinate

	for m, mol in mols:
		for t, type in mol.types:
			molecule_name.append(mol.molid)
			atom_index.append(t)
			atom.append(p_table(type))
			x.append(mol.xyz[t][0])
			y.append(mol.xyz[t][1])
			z.append(mol.xyz[t][2])

	atoms = {	'molecule_name': molecule_name,
				'atom_index': atom_index,
				'atom': atom,
				'x': x,
				'y': y,
				'z': z
			}

	structure_dict = BCAI.make_strucure_dict(atoms)

	BCAI.enhance_structure_dict(structure_dict)

	BCAI.enhance_atoms(atoms, structure_dict)

	# construct dataframe as BCAI requires from mols
	# atoms has: molecule_name, atom, labeled atom,
	id = []				# number
	molecule_name = [] 	# molecule name
	atom_index_0 = []	# atom index for atom 1
	atom_index_1 = []	# atom index for atom 2
	type = []			# coupling type
	scalar_coupling_constant = []	# coupling value

	i = -1
	for m, mol in mols:
		for t, type in mol.types:
			for t2, type2 in mol.types:
				if t == t2:
					continue
				if type != target[1] and type2 != target[2]:
					continue
				if mol.coupling_len[t][t2] != target[0]:
					continue

				i += 1
				id.append(i)
				molecule_name.append(mol.molid)
				atom_index_0.append(t)
				atom_index_1.append(t2)
				type.append(targetflag)
				scalar_coupling_constant.append(mol.coupling[t][t2])

	bonds = {	'id': id,
				'molecule_name': molecule_name,
				'atom_index_0': atom_index,
				'atom_index_1': atom,
				'type': type,
				'scalar_coupling_constant': coupling
			}


	BCAI.enhance_bonds(bonds, structure_dict)
	bonds = BCAI.add_all_pairs(bonds, structure_dict) # maybe replace this
	triplets = BCAI.make_triplets(bonds["molecule_name"].unique(), structure_dict)
	quadruplets = make_quadruplets(bonds["molecule_name"].unique(),structure_dict)

	atoms = pd.DataFrame(atoms)
	bonds = pd.DataFrame(bonds)
	triplets = pd.DataFrame(triplets)
	qudadruplets = pd.DataFrame(quadruplets)

	atoms.sort_values(['molecule_name','atom_index'],inplace=True)
	bonds.sort_values(['molecule_name','atom_index_0','atom_index_1'],inplace=True)
	triplets.sort_values(['molecule_name','atom_index_0','atom_index_1','atom_index_2'],inplace=True)

	embeddings = BCAI.add_embedding(atoms, bonds, triplets, quadruplets)
	means, stds = get_scaling(bonds)
	bonds = add_scaling(bonds, means, stds)


	x = BCAI.create_dataset(atoms, bonds, triplets, quads, labeled = True, max_count = 10**10)
	y = x[-1]

	print(x)
	sys.exit(0)

	return x, y, r










































###

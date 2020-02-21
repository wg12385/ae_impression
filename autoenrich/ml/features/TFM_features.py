# Copyright 2020 Will Gerrard
#This file is part of autoenrich.

#autoenrich is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoenrich is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.

# get features for TransForMer models
from autoenrich.reference.periodic_table import Get_periodic_table
from autoenrich.util.flag_handler.hdl_targetflag import flag_to_target
from autoenrich.molecule.nmrmol import nmrmol
from autoenrich.util.file_gettype import get_type

import autoenrich.ml.features.BCAI_calc.mol_graph_setup as BCAI

import numpy as np
import pandas as pd
import sys
import pickle
import gzip
from tqdm import tqdm

# make BCAI features
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
	conns = []

	mol_order = []
	m = -1
	for molrf in tqdm(mols, desc='Constructing atom dictionary'):
		m += 1
		if len(mols) > 2000:
			mol = nmrmol(molid=molrf[1])

			if molrf[2] == '':
				ftype = get_type(molrf[2])
			else:
				ftype = molrf[2]
			mol.read_nmr(molrf[0], ftype)
		else:
			mol = molrf
		mol_order.append(mol.molid)
		for t, type in enumerate(mol.types):
			molecule_name.append(mol.molid)
			atom_index.append(t)
			atom.append(p_table[type])
			x.append(mol.xyz[t][0])
			y.append(mol.xyz[t][1])
			z.append(mol.xyz[t][2])
			conns.append(mol.conn[t])

	atoms = {	'molecule_name': molecule_name,
				'atom_index': atom_index,
				'atom': atom,
				'x': x,
				'y': y,
				'z': z,
				'conn': conns,
			}

	atoms = pd.DataFrame(atoms)
	structure_dict = BCAI.make_structure_dict(atoms)

	BCAI.enhance_structure_dict(structure_dict)

	BCAI.enhance_atoms(atoms, structure_dict)

	# construct dataframe as BCAI requires from mols
	# atoms has: molecule_name, atom, labeled atom,
	id = []				# number
	molecule_name = [] 	# molecule name
	atom_index_0 = []	# atom index for atom 1
	atom_index_1 = []	# atom index for atom 2
	cpltype = []			# coupling type
	coupling = []	# coupling value
	r = []
	y = []

	i = -1
	m = -1
	for molrf in tqdm(mols, desc='Constructing bond dictionary'):
		m += 1
		if len(mols) > 2000:
			mol = nmrmol(molid=molrf[1])

			if molrf[2] == '':
				ftype = get_type(molrf[2])
			else:
				ftype = molrf[2]
			mol.read_nmr(molrf[0], ftype)
		else:
			mol = molrf

		for t, type in enumerate(mol.types):
			for t2, type2 in enumerate(mol.types):
				if t == t2:
					continue
				if not ( type == target[1] and type2 == target[2] ):
					continue
				if mol.coupling_len[t][t2] != target[0]:
					continue

				i += 1
				id.append(i)
				molecule_name.append(mol.molid)
				atom_index_0.append(t)
				atom_index_1.append(t2)

				TFM_flag = targetflag[2] + '-' + targetflag[3] + '_' + targetflag[0] + '.0'

				cpltype.append(TFM_flag)

				coupling.append(mol.coupling[t][t2])

				if np.isnan(mol.coupling[t][t2]):
					print(mol.molid)

				y.append(mol.coupling[t][t2])
				r.append([mol.molid, t, t2])

	bonds = {	'id': id,
				'molecule_name': molecule_name,
				'atom_index_0': atom_index_0,
				'atom_index_1': atom_index_1,
				'type': cpltype,
				'scalar_coupling_constant': coupling
			}

	#print(len(id), len(molecule_name), len(atom_index), len(atom))

	bonds = pd.DataFrame(bonds)

	bonds = BCAI.enhance_bonds(bonds, structure_dict)

	bonds = BCAI.add_all_pairs(bonds, structure_dict) # maybe replace this
	triplets = BCAI.make_triplets(bonds["molecule_name"].unique(), structure_dict)

	atoms = pd.DataFrame(atoms)
	bonds = pd.DataFrame(bonds)
	triplets = pd.DataFrame(triplets)

	atoms.sort_values(['molecule_name','atom_index'],inplace=True)
	bonds.sort_values(['molecule_name','atom_index_0','atom_index_1'],inplace=True)
	triplets.sort_values(['molecule_name','atom_index_0','atom_index_1','atom_index_2'],inplace=True)

	embeddings, atoms, bonds, triplets = BCAI.add_embedding(atoms, bonds, triplets)
	bonds.dropna()
	atoms.dropna()
	means, stds = BCAI.get_scaling(bonds)
	bonds = BCAI.add_scaling(bonds, means, stds)

	x = BCAI.create_dataset(atoms, bonds, triplets, labeled = True, max_count = 10**10, mol_order=mol_order)

	return x, y, r, mol_order










































###

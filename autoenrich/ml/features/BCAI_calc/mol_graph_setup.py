#!/usr/bin/env python

## Copyright (c) 2017 Robert Bosch GmbH
## All rights reserved.
##
## This source code is licensed under the MIT license found in the
## LICENSE file in the root directory of this source tree.

# Edited, Will Gerrard, 2020, for use in autoenrich

import collections
import gzip
import itertools
import json
import os
import pickle
import sys

import numpy as np
import pandas as pd
import rdkit
import autoenrich.ml.features.BCAI_calc.xyz2mol as x2m
from autoenrich.reference.periodic_table import Get_periodic_table
from tqdm import tqdm

# Due to some compatibility issues between rdkit/pybel and torch, we have to load them as needed.
# Rules are meant to be broken, including best-programming practices :)


bond_order_dict = { rdkit.Chem.rdchem.BondType.SINGLE: 1,
					rdkit.Chem.rdchem.BondType.AROMATIC: 1.5,
					rdkit.Chem.rdchem.BondType.DOUBLE: 2,
					rdkit.Chem.rdchem.BondType.TRIPLE: 3}

atomic_num_dict = { 'H':1, 'C':6, 'N':7, 'O':8, 'F':9 }
# These were mistaken or too small datasets, so we are relabeling them.
classification_corrections = {
					  '1JHN_2_2_1_1':'1JHN_3_2_2_1',
					  '3JHN_4.5_3_1.5_1.5':'3JHN_4_3_1.5_1.5',
					  '2JHC_3_3_1_1':'2JHC_4_3_2_1',
					  '3JHC_3_3_1_1':'3JHC_4_3_2_1',
					  '3JHC_4_2_2_2':'3JHC_4_2_3_1'}
# These have less than 1000 between train and test, so we will drop the subtypes
small_longtypes = {'2JHN_4.5_2_3_1.5', '3JHN_4_2_3_1', '2JHN_4_2_3_1',
				   '2JHN_4.5_3_1.5_1.5', '2JHN_4_3_2_1', '3JHN_4_4_1_1',
				   '3JHN_4_3_2_1', '2JHN_4_4_1_1', '3JHN_4.5_2_3_1.5',
				   '2JHN_4_2_2_2', '3JHN_4_2_2_2', '1JHN_4_3_2_1',
				   '1JHN_4_4_1_1', '2JHN_3_1_3_0'}
(MAX_ATOM_COUNT,MAX_BOND_COUNT,MAX_TRIPLET_COUNT,MAX_QUAD_COUNT) = (29, 406, 54, 117)

def make_atom_df(mols):

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

	return atoms

def make_bonds_df(mols):

		p_table = Get_periodic_table()

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

			moly = []
			molr = []

			for t, type in enumerate(mol.types):
				for t2, type2 in enumerate(mol.types):
					if t == t2:
						continue

					TFM_flag = str(mol.coupling_len[t][t2]) + 'J' + p_table[type] + p_table[type2]

					#if TFM_flag != flag and flag != 'all':
					#	continue

					i += 1
					id.append(i)
					molecule_name.append(mol.molid)
					atom_index_0.append(t)
					atom_index_1.append(t2)


					cpltype.append(TFM_flag)

					coupling.append(mol.coupling[t][t2])

					moly.append(mol.coupling[t][t2])
					molr.append([mol.molid, t, t2])

			y.append(moly)
			r.append(molr)

		bonds = {	'id': id,
					'molecule_name': molecule_name,
					'atom_index_0': atom_index_0,
					'atom_index_1': atom_index_1,
					'type': cpltype,
					'scalar_coupling_constant': coupling
				}

		bonds = pd.DataFrame(bonds)

		return bonds

def make_structure_dict(atoms_dataframe):
	"""Convert from structures.csv output to a dictionary data storage.

	Args:
		atoms_dataframe: The dataframe corresponding to structures.csv

	Returns:
		dict: Mapping of molecule name to molecule properties.

	"""
	atoms = atoms_dataframe.sort_values(["molecule_name", "atom_index"])  # ensure ordering is consistent
	# Make a molecule-based dictionary of the information
	structure_dict = collections.defaultdict(lambda: {"symbols":[],"positions":[],"conn":[]})
	for index,row in atoms.iterrows():
		structure_dict[row["molecule_name"]]["symbols"].append(row["atom"])
		structure_dict[row["molecule_name"]]["positions"].append([row["x"],row["y"],row["z"]])
		structure_dict[row["molecule_name"]]["conn"].append(row["conn"])

	return structure_dict


def enhance_structure_dict(structure_dict):
	"""Add derived information to the structure dictionary.

	Args:
		structure_dict: Output of :func:`make_structure_dict`.

	Returns:
		dict: The same, modified in-place, with derived information (e.g. atom distances).

	Caution: If torch is imported at the same time as this is run, you may get a segmentation fault. Complain to pybel or rdkit, I suppose.
	"""

	import pybel

	atomic_num_dict = { 'H':1, 'C':6, 'N':7, 'O':8, 'F':9 }

	for molecule_name in tqdm(structure_dict, desc='Enhancing atom dictionary'):

		# positions - array (N,3) of Cartesian positions
		molecule = structure_dict[molecule_name]
		positions = np.array(molecule['positions'])
		conn = np.array(molecule['conn'])
		n_atom = positions.shape[0]
		molecule['positions'] = positions

		# distances - array (N,N) of distances between atoms
		pos1 = np.tile(positions, (n_atom,1,1) )
		pos2 = np.transpose(pos1, (1,0,2) )
		dist = np.linalg.norm(pos1 - pos2, axis=-1)
		molecule['distances'] = dist

		# angle - array (N,) of angles to the 2 closest atoms
		sorted_j = np.argsort(dist, axis=-1)
		relpos1 = positions[sorted_j[:,1],:] - positions[sorted_j[:,0],:]
		relpos2 = positions[sorted_j[:,2],:] - positions[sorted_j[:,0],:]
		cos = np.sum(relpos1*relpos2,axis=1) / (np.linalg.norm(relpos1,axis=1) * np.linalg.norm(relpos2,axis=1))
		angle = np.arccos( np.clip(cos,-1.0,1.0) ).reshape((n_atom,1)) / np.pi
		molecule['angle'] = angle[:,0]

		# bond orders - array (N,N) of the bond order (0 for no chemical bond)
		molecule['bond_orders'] = np.zeros((n_atom,n_atom))
		atomicNumList = [atomic_num_dict[symbol] for symbol in molecule['symbols']]

		for atom0 in range(len(molecule['symbols'])):
			for atom1 in range(len(molecule['symbols'])):
				try:
					molecule['bond_orders'][atom0,atom1] = conn[atom0][atom1]
					molecule['bond_orders'][atom1,atom0] = conn[atom1][atom0]
				except Exception as e:
					print(e)
					print(molecule_name)
					print(atom0, atom1)
					print('positions:', molecule['positions'].shape, positions.shape)
					print('conn:', len(molecule['conn']), len(molecule['conn'][0]), conn.shape)
					print('conn[0]:', conn[0].shape)
					print('bond_orders:', molecule['bond_orders'].shape)
					print(n_atom)
					sys.exit(0)

		# Supplementary information for tagging:
		# top_bonds: (N,4 or less) bond orders of the top 4 bonds, for each atom
		# bond_ids: (N,4): Label the atom with the following 4 linear transform of top_bonds:
		#   * total num bonds (valence), counting double as 2
		#   * total num bonded neighbors, counting double as 1
		#   * largest order
		#   * second largest order.
		molecule['top_bonds'] = np.sort(molecule['bond_orders'],axis=-1)[:,-1:-5:-1]
		molecule['bond_ids'] = np.hstack((molecule['top_bonds'].sum(axis=-1)[:,np.newaxis],
										  np.sum(molecule['top_bonds']>1e-3,axis=-1)[:,np.newaxis],
										  molecule['top_bonds'][:,:2]))
		# long_symbols (N,) string relabel of the symbol straight from bond_ids
		molecule['long_symbols'] = ['_'.join([
			molecule['symbols'][i]]+[str(x) for x in molecule['bond_ids'][i]])
									for i in range(n_atom)]
		chem_bond_atoms = [sorted([molecule['symbols'][i] for i in molecule['bond_orders'][atom_index].nonzero()[0]])
						   for atom_index in range(n_atom)]
		molecule['sublabel_atom'] = ['-'.join([molecule['long_symbols'][atom_index]]+chem_bond_atoms[atom_index])
									for atom_index in range(n_atom)]

		# pybel information. I think we only end up using Gastiger charges.
		# Each of these is (N,) arrays
		# Convert to xyz string for pybel's I/O
		xyz = str(n_atom)+'\n\n' + '\n'.join([ ' '.join( [
				str(molecule['symbols'][i]),
				str(molecule['positions'][i,0]),
				str(molecule['positions'][i,1]),
				str(molecule['positions'][i,2])] )
				for i in range(n_atom)])

		mol = pybel.readstring('xyz',xyz)
		molecule['charges'] = [mol.atoms[i].partialcharge for i in range(n_atom)]

	return structure_dict


def enhance_atoms(atoms_dataframe,structure_dict):
	"""Enhance the atoms dataframe by including derived information.

	Args:
		atoms_dataframe: Pandas dataframe read from structures.csv.
		structure_dict: Output of :func:`make_structure_dict`, after running :func:`enhance_structure_dict`.

	Returns:
		pandas.DataFrame: Same dataframe, modified in-place, with derived information added.

	"""
	assert int(atoms_dataframe.groupby("molecule_name").count().max()[0]) <= MAX_ATOM_COUNT
	for key in tqdm(['distances','angle', 'bond_orders', 'top_bonds', 'bond_ids', 'long_symbols','sublabel_atom',
				'charges'], desc='Enhancing atoms'):
		newkey = key if key[-1]!='s' else key[:-1]
		atoms_dataframe[newkey] = atoms_dataframe.apply(lambda x:
														structure_dict[x['molecule_name']][key][x['atom_index']],
														axis=1)
		atoms_dataframe.rename(columns={'long_symbol':'labeled_atom'},inplace=True)
	return atoms_dataframe


def enhance_bonds(bond_dataframe, structure_dict, flag='3JHH'):
	"""Enhance the bonds dataframe by including derived information.

	Args:
		bond_dataframe: Pandas dataframe read from train.csv or test.csv.
		structure_dict: Output of :func:`make_structure_dict`, after running :func:`enhance_structure_dict`.

	Returns:
		pandas.DataFrame: Same dataframe, modified in-place, with derived information added.

	"""
	bond_dataframe.sort_values(['molecule_name','atom_index_0','atom_index_1'],inplace=True)
	assert int(bond_dataframe.groupby("molecule_name").count().max()[0]) <= MAX_BOND_COUNT
	new_columns = collections.defaultdict(list)
	for index,row in bond_dataframe.iterrows():
		molecule_name, iatom0, iatom1 = row['molecule_name'],row['atom_index_0'],row['atom_index_1']
		if 'predict' not in structure_dict[molecule_name]:
			structure_dict[molecule_name]['predict'] = structure_dict[molecule_name]['bond_orders'] * 0
		structure_dict[molecule_name]['predict'][iatom0,iatom1] = 1
		structure_dict[molecule_name]['predict'][iatom1,iatom0] = 1
		long_symbols = [structure_dict[molecule_name]['long_symbols'][x] for x in [iatom0,iatom1]]

		# labeled_type
		if all([x[0]=='H' for x in long_symbols]):
			lt = row['type']
		else:
			ls = [x for x in long_symbols if x[0]!='H'][0]
			lt = row["type"] + ls[1:].replace('.0','')
			if lt in classification_corrections:
				lt = classification_corrections[lt]
			if lt in small_longtypes:
				lt = lt.split('_')[0]

		new_columns["labeled_type"].append(lt)

		# sublabeled type
		new_columns["sublabel_type"].append(row['type'] + '-'+ '-'.join(sorted(long_symbols)))
		# bond order
		new_columns["bond_order"].append(structure_dict[molecule_name]['bond_orders'][iatom0,iatom1])

		if lt == flag:
			new_columns["predict"].append(1)
		else:
			new_columns["predict"].append(0)
	for key in new_columns:
		bond_dataframe[key] = new_columns[key]
	return bond_dataframe

def make_triplets(molecule_list,structure_dict):
	"""Make the triplet dataframe.

	Args:
		molecule_list: List of molecules to generate.
		structure_dict: Output of :func:`make_structure_dict`, after running :func:`enhance_structure_dict`.

	Returns:
		pandas.DataFrame: New dataframe, with triplets and related information. The convention is the bond looks like 1-0-2, where 0 is the central atom.

	"""
	new_data = collections.defaultdict(list)
	for molecule_name in molecule_list:
		molecule = structure_dict[molecule_name]
		bond_orders = molecule['bond_orders']
		short = molecule['symbols']
		long = molecule['long_symbols']
		for i, atom_bond_order in enumerate(bond_orders):
			connection_indices = atom_bond_order.nonzero()[0]
			pairs = itertools.combinations(connection_indices,2)
			for pair in pairs:
				j, k = pair[0], pair[1]

				atom0_short = short[i] + long[i].split('_')[2]
				atom1_short = short[j] + long[j].split('_')[2]
				atom2_short = short[k] + long[k].split('_')[2]
				atom0_long = long[i]
				atom1_long = long[j]
				atom2_long = long[k]
				#labels = ['-'.join([atom1_short,str(atom_bond_order[j])]),
				#          '-'.join([atom2_short,str(atom_bond_order[k])])]
				labels = [atom1_short,atom2_short]
				labels.sort()
				label = '-'.join([atom0_short]+labels)
				#sublabels = ['-'.join([atom1_long,str(atom_bond_order[j])]),
				#             '-'.join([atom2_long,str(atom_bond_order[k])])]
				sublabels = [atom1_long,atom2_long]
				sublabels.sort()
				sublabel = '-'.join([atom0_long]+sublabels)
				r10 = molecule['positions'][j] - molecule['positions'][i]
				r20 = molecule['positions'][k] - molecule['positions'][i]
				angle = np.sum(r10*r20) / (np.linalg.norm(r10)*np.linalg.norm(r20))
				angle = np.arccos( np.clip(angle,-1.0,1.0) )
				row = {'molecule_name':molecule_name,'atom_index_0':i,'atom_index_1':j,'atom_index_2':k,
					  'label':label,'sublabel':sublabel,'angle':angle}
				for k,v in row.items():
					new_data[k].append(v)
	ans = pd.DataFrame(new_data)
	ans.sort_values(['molecule_name','atom_index_0','atom_index_1','atom_index_2'])
	assert int(ans.groupby("molecule_name").count().max()[0]) <= MAX_TRIPLET_COUNT
	return ans

def write_csv(directory,label,atoms,bonds,triplets):
	"""Write the relevant dataframes to a CSV file.

	Args:
		directory: Directory to write to.
		label (str): How to label the files, e.g. test or train.
		atoms: Pandas dataframe read from structures.csv, after running :func:`enhance_atoms`.
		bonds: Pandas dataframe read from train.csv or test.csv, after running :func:`enhance_bonds`.
		triplets: Pandas dataframe created by :func:`make_triplets`.

	Returns:
		None

	"""
	filename = os.path.join(directory,'new_big_{}.csv.bz2')
	if atoms is not None and len(atoms):
		atoms = atoms.sort_values(["molecule_name",'atom_index'])
		for i in range(4):
			atoms["top_bond_{}".format(i)] = [x[i] if len(x)>i else 0.0 for x in atoms["top_bond"].values]
		for i in ["x","y","z"]:
			atoms[i] = atoms[i].values.round(10)
		renames = {k:k[:-1] for k in atoms.columns if k[-1]=='s'}
		renames.update({'long_symbols':'labeled_atom'})
		atoms = atoms.rename(columns=renames)
		atoms.to_csv(filename.format('structures'),index=False,columns=
			'molecule_name,atom_index,atom,x,y,z,labeled_atom,angle,top_bond_0,top_bond_1,top_bond_2,top_bond_3,sublabel_atom,charge,spin,heavyvalence,heterovalence,valence,hyb_type'.split(','))
	if bonds is not None and len(bonds):
		bonds = bonds.reset_index()
		bond_columns = 'id,molecule_name,atom_index_0,atom_index_1,type,scalar_coupling_constant,labeled_type,sublabel_type,bond_order,predict'.split(',')
		if 'scalar_coupling_constant' not in bonds.columns:
			bond_columns = [x for x in bond_columns if x!='scalar_coupling_constant']
		bonds = bonds.sort_values(["predict","molecule_name",'atom_index_0','atom_index_1'],
								  ascending=[False,True,True,True])
		bonds.to_csv(filename.format(label),index=False,columns=bond_columns)
	if triplets is not None and len(triplets):
		triplets = triplets.sort_values(["molecule_name",'atom_index_0','atom_index_1','atom_index_2'])
		triplets.to_csv(filename.format(label+'_triplets'),index=False,columns=
			'molecule_name,atom_index_0,atom_index_1,atom_index_2,label,sublabel,angle'.split(','))

def _create_embedding(series):
	"""Create a one-hot encoding embedding.

	Args:
		series: A DataFrame series (column).

	Returns:
		dict: Mapping of the entries (or "<None>") to the index number.

	"""
	types = sorted(series.unique().tolist())
	assert "<None>" not in types
	emb_index = dict(zip(["<None>"] + types , range(len(types)+1)))
	return emb_index


def add_embedding(atoms,bonds,triplets,embeddings=None):
	"""Add embedding indices to the dataframes.

	Args:
		atoms: Pandas dataframe read from structures.csv, after running :func:`enhance_atoms`.
		bonds: Pandas dataframe read from train.csv or test.csv, after running :func:`enhance_bonds`.
		triplets: Pandas dataframe created by :func:`make_triplets`.
		embeddings (dict or None): If None, we create a new embedding (e.g. train data), otherwise we use the given embeddings thar are output by :func:`add_embedding` (e.g. test data).

	Returns:
		dict: The embedding dictionary that can be passed to this function for using the same embedding on a new dataset.

	"""
	# Add the embedding info to the dataframes.
	atoms["type_0"] = atoms["atom"]
	atoms["type_1"] = atoms["labeled_atom"].apply(lambda x : x[:5])
	atoms["type_2"] = atoms["labeled_atom"]
	bonds["type_0"] = bonds["type"]
	bonds["type_1"] = bonds["labeled_type"]
	bonds["type_2"] = bonds["sublabel_type"]
	triplets["type_0"] = triplets["label"].apply(lambda x : x[0] + x[5] + x[10])
	triplets["type_1"] = triplets["label"]
	if embeddings is None:
		embeddings = {}
		embeddings.update({('atom',t):_create_embedding(atoms["type_" + str(t)]) for t in range(3)})
		embeddings.update({('bond',t):_create_embedding(bonds["type_" + str(t)]) for t in range(3)})
		embeddings.update({('triplet',t):_create_embedding(triplets["type_" + str(t)]) for t in range(2)})
	for t in range(3):
		atoms["type_index_" + str(t)] = atoms["type_" + str(t)].apply(lambda x : embeddings[('atom',t)][x])
	for t in range(3):
		bonds["type_index_" + str(t)] = bonds["type_" + str(t)].apply(lambda x : embeddings[('bond',t)][x])
	for t in range(2):
		triplets["type_index_" + str(t)] = triplets["type_" + str(t)].apply(lambda x : embeddings[('triplet',t)][x])

	return embeddings, atoms, bonds, triplets


def get_scaling(bonds_train):
	"""Get the mean/std scaling factors for each ``labeled_type``.

	Args:
		bonds_train: The training data that we can use to set the values.

	Returns:
		tuple: Mean and std dicts, mapping labeled_type to scalar_coupling_constant mean/std.

	"""
	# Get the mean/std scaling factors
	means = bonds_train.groupby("labeled_type").mean()["scalar_coupling_constant"].to_dict()
	stds = bonds_train.groupby("labeled_type").std()["scalar_coupling_constant"].to_dict()
	# stds of 0 were causing NaNs in training, assuming this workflow originally assumed at least one of
	# every coupling type would be present
	for key in stds:
		if np.isnan(stds[key]) or stds[key] == 0:
			stds[key] = 1

	return means,stds


def add_scaling(bonds,means,stds):
	"""Add the scaling information to the bonds dataframe.

	Args:
		bonds (pd.DataFrame): The dataframe of the bonds, after :func:`enhance_bonds`.
		means (dict): Output of :func:`get_scaling`.
		stds (dict): Output of :func:`get_scaling`.

	Returns:
		pd.DataFrame: Same dataframe, with added columns.

	"""
	# Add mean/std scaling factors to bonds dataframe
	bonds["sc_mean"] = bonds["labeled_type"].apply(lambda x : means[x])
	bonds["sc_std"] = bonds["labeled_type"].apply(lambda x : stds[x])
	if "scalar_coupling_constant" in bonds.columns:
		bonds["sc_scaled"] = (bonds["scalar_coupling_constant"] - bonds["sc_mean"]) / bonds["sc_std"]
	return bonds


def create_dataset(atoms, bonds, triplets, labeled = True, max_count = 10**10, mol_order=[]):
	"""Create the python loaders, which we can pkl to a file for batching.

	Args:
		atoms: Pandas dataframe read from structures.csv, after running :func:`enhance_atoms`.
		bonds: Pandas dataframe read from train.csv or test.csv, after running :func:`enhance_bonds`.
		triplets: Pandas dataframe created by :func:`make_triplets`.
		labeled (bool): Whether this is train data, labeled with the y value.
		max_count (int): Maximum number of entries; useful for testing.

	Returns:
		tuple: With the following entries

			* x_index: (M,) Index of the molecule.
			* x_atom: (M,N,3) Atom type index.
			* x_atom_pos: (M,N,5) Atom position (3), closest-atom angle (1), and partial charge (1).
			* x_bond: (M,B,5) Bond type index (3), Atom index (2) corresponding to the bond.
			* x_bond_dist: (M,B) Distance of the bond.
			* x_triplet: (N,P,7): Triplet type (2), Atom index (3), Bond index (2) corresponding to the triplet.
			* x_triplet_angle: (N,P) Triplet angle.
			* y_bond_scalar_coupling: (N,M,4) of the scalar coupling constant, type mean, type std, and whether it should be predicted.

	"""
	import torch
	from tqdm import tqdm
	# create mapping from molecule names to indices
	mol_unique = sorted(bonds["molecule_name"].unique().tolist())
	index = dict(zip(mol_unique, range(len(mol_unique))))
	atoms = atoms.set_index("molecule_name")
	bonds = bonds.set_index("molecule_name")
	triplets = triplets.set_index("molecule_name")

	max_count = M = min(max_count, len(index))

	x_index = torch.arange(M, dtype=torch.long)
	x_atom = torch.zeros(M, MAX_ATOM_COUNT, 3, dtype=torch.long)
	x_atom_pos = torch.zeros(M, MAX_ATOM_COUNT, 5)
	x_bond = torch.zeros(M, MAX_BOND_COUNT, 5, dtype=torch.long)
	x_bond_dist = torch.zeros(M, MAX_BOND_COUNT)
	x_triplet = torch.zeros(M, MAX_TRIPLET_COUNT, 7, dtype=torch.long)
	x_triplet_angle = torch.zeros(M, MAX_TRIPLET_COUNT)

	y_bond_scalar_coupling = torch.zeros(M, MAX_BOND_COUNT, 4)

	i = -1
	for k in tqdm(mol_order, desc='Constructing training dataset'):
		i += 1
		if i >= M:
			break
		mol_atoms = atoms.loc[[k]]
		mol_bonds = bonds.loc[[k]]
		mol_real_bonds = mol_bonds[(mol_bonds["predict"]==1) | (mol_bonds["bond_order"]>0)]
		mol_fake_bonds = mol_bonds[(mol_bonds["predict"]==0) & (mol_bonds["bond_order"]==0)]
		mol_triplets = triplets.loc[[k]]

		n = mol_atoms.shape[0]
		m = mol_bonds.shape[0]
		mr = mol_real_bonds.shape[0]
		mf = mol_fake_bonds.shape[0]
		p = mol_triplets.shape[0]
		assert mr + mf == m, "Real + fake bonds != number of bonds?"
		assert mr < MAX_BOND_COUNT, "The number of real bonds is SMALLER than the MAX_BOND_COUNT"

		# STEP 1: Atoms
		for t in range(3):
			x_atom[i,:n,t] = torch.tensor(mol_atoms["type_index_" + str(t)].values)
		x_atom_pos[i,:n,:3] = torch.tensor(mol_atoms[["x", "y", "z"]].values)
		x_atom_pos[i,:n,3] = torch.tensor(mol_atoms["angle"].values)
		x_atom_pos[i,:n,4] = torch.tensor(mol_atoms["charge"].values)

		# STEP 2: Real bonds
		for t in range(3):
			x_bond[i,:mr,t] = torch.tensor(mol_real_bonds["type_index_" + str(t)].values)
		x_bond[i,:mr,3] = torch.tensor(mol_real_bonds["atom_index_0"].values)
		x_bond[i,:mr,4] = torch.tensor(mol_real_bonds["atom_index_1"].values)

		idx1 = torch.tensor(mol_real_bonds["atom_index_0"].values)
		idx2 = torch.tensor(mol_real_bonds["atom_index_1"].values)
		x_bond_dist[i,:mr] = ((x_atom_pos[i,idx1,:3] - x_atom_pos[i,idx2,:3])**2).sum(1)

		if mf > 0:
			# STEP 3: Fake bonds
			fidx1 = torch.tensor(mol_fake_bonds["atom_index_0"].values)
			fidx2 = torch.tensor(mol_fake_bonds["atom_index_1"].values)
			fdists = ((x_atom_pos[i,fidx1,:3] - x_atom_pos[i,fidx2,:3])**2).sum(1)   # Length mf
			argsort_fdists = torch.argsort(fdists)
			top_count = min(MAX_BOND_COUNT - mr, mf)

			for t in range(3):
				x_bond[i,mr:mr+top_count,t] = torch.tensor(mol_fake_bonds["type_index_" + str(t)].values)[argsort_fdists][:top_count]
			x_bond[i,mr:mr+top_count,3] = torch.tensor(mol_fake_bonds["atom_index_0"].values)[argsort_fdists][:top_count]
			x_bond[i,mr:mr+top_count,4] = torch.tensor(mol_fake_bonds["atom_index_1"].values)[argsort_fdists][:top_count]
			x_bond_dist[i,mr:mr+top_count] = fdists[argsort_fdists][:top_count]

		# STEP 4: Triplets
		for t in range(2):
			x_triplet[i,:p,t] = torch.tensor(mol_triplets["type_index_" + str(t)].values)
		x_triplet[i,:p,2] = torch.tensor(mol_triplets["atom_index_0"].values)
		x_triplet[i,:p,3] = torch.tensor(mol_triplets["atom_index_1"].values)
		x_triplet[i,:p,4] = torch.tensor(mol_triplets["atom_index_2"].values)

		x_triplet_angle[i,:p] = torch.tensor(mol_triplets["angle"].values)
		lookup = dict(zip(mol_real_bonds["atom_index_0"].apply(str) + "_" + mol_real_bonds["atom_index_1"].apply(str),
						  range(mol_real_bonds.shape[0])))
		lookup.update(dict(zip(mol_real_bonds["atom_index_1"].apply(str) + "_" + mol_real_bonds["atom_index_0"].apply(str),
						   range(mol_real_bonds.shape[0]))))

		b_idx1 = (mol_triplets["atom_index_0"].apply(str) + "_" +
				  mol_triplets["atom_index_1"].apply(str)).apply(lambda x : lookup[x])
		b_idx2 = (mol_triplets["atom_index_0"].apply(str) + "_" +
				  mol_triplets["atom_index_2"].apply(str)).apply(lambda x : lookup[x])

		x_triplet[i,:p,5] = torch.tensor(b_idx1.values)
		x_triplet[i,:p,5] = torch.tensor(b_idx2.values)

		if labeled:
			y_bond_scalar_coupling[i,:mr, 0] = torch.tensor(mol_real_bonds["scalar_coupling_constant"].values)
		else:
			y_bond_scalar_coupling[i,:mr, 0] = torch.tensor(mol_real_bonds["id"].values)
		y_bond_scalar_coupling[i,:mr, 1] = torch.tensor(mol_real_bonds["sc_mean"].values)
		y_bond_scalar_coupling[i,:mr, 2] = torch.tensor(mol_real_bonds["sc_std"].values)
		y_bond_scalar_coupling[i,:mr, 3] = torch.tensor(mol_real_bonds["predict"].values).float()  # binary tensor (1s to be predicted)

	return x_index, x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, y_bond_scalar_coupling

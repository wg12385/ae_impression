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

import tracemalloc

import numpy as np
import pandas as pd
import sys
import pickle
import gzip
from tqdm import tqdm

import torch

def get_size(obj, seen=None):
	"""Recursively finds size of objects"""
	size = sys.getsizeof(obj)
	if seen is None:
		seen = set()
	obj_id = id(obj)
	if obj_id in seen:
		return 0
	# Important mark as seen *before* entering recursion to gracefully handle
	# self-referential objects
	seen.add(obj_id)
	if isinstance(obj, dict):
		size += sum([get_size(v, seen) for v in obj.values()])
		size += sum([get_size(k, seen) for k in obj.keys()])
	elif hasattr(obj, '__dict__'):
		size += get_size(obj.__dict__, seen)
	elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
		size += sum([get_size(i, seen) for i in obj])
	return size

def save_dataset(Dset, y, r, mol_order):
	p = np.random.permutation(Dset[0].shape[0])
	idx_train = torch.cat([torch.tensor(p[:int(0.6*len(p))]), torch.tensor(p[int(0.8*len(p)):])])
	idx_val = torch.tensor(p[int(0.6*len(p)):int(0.8*len(p))])

	D_train = tuple([d[idx_train] for d in Dset])
	D_val = tuple([d[idx_val] for d in Dset])

	train_file = "training_data/dataset_features_x.pkl.gz"
	with gzip.open(train_file, "wb") as f:
		pickle.dump(Dset, f, protocol=4)
	x = [train_file, train_file]

	train_file = "training_data/dataset_features_r.pkl.gz"
	with gzip.open(train_file, "wb") as f:
		pickle.dump(r, f, protocol=4)
	r = [train_file, train_file]

	train_file = "training_data/dataset_features_y.pkl.gz"
	with gzip.open(train_file, "wb") as f:
		pickle.dump(y, f, protocol=4)
	y = [train_file, train_file]

	train_file = "training_data/dataset_features_m.pkl.gz"
	with gzip.open(train_file, "wb") as f:
		pickle.dump(mol_order, f, protocol=4)
	mol_order = [train_file, train_file]

	return x, y, r, mol_order


def save_split_dataset(Dset, y, r, mol_order):

	snapshot = tracemalloc.take_snapshot()
	top_stats = snapshot.statistics('lineno')
	print("12[ Top 2 ]")
	for stat in top_stats[:2]:
		print(stat)


	p = np.random.permutation(Dset[0].shape[0])
	idx_train = torch.cat([torch.tensor(p[:int(0.6*len(p))]), torch.tensor(p[int(0.8*len(p)):])])
	idx_val = torch.tensor(p[int(0.6*len(p)):int(0.8*len(p))])

	D_train = tuple([d[idx_train] for d in Dset])
	D_val = tuple([d[idx_val] for d in Dset])


	train_file = "training_data/dataset_features_xtrain.pkl.gz"
	val_file = "training_data/dataset_features_xval.pkl.gz"
	with gzip.open(train_file, "wb") as f:
		pickle.dump(D_train, f, protocol=4)
	with gzip.open(val_file, "wb") as f:
		pickle.dump(D_val, f, protocol=4)
	x = [train_file, val_file]

	idx_train = [np.asscalar(id.numpy()) for id in idx_train]
	idx_val = [np.asscalar(id.numpy()) for id in idx_val]

	r_train = [r[id] for id in idx_train]
	r_val = [r[id] for id in idx_val]
	train_file = "training_data/dataset_features_rtrain.pkl.gz"
	val_file = "training_data/dataset_features_rval.pkl.gz"
	with gzip.open(train_file, "wb") as f:
		pickle.dump(r_train, f, protocol=4)
	with gzip.open(val_file, "wb") as f:
		pickle.dump(r_val, f, protocol=4)
	r = [train_file, val_file]


	y_train = [y[id] for id in idx_train]
	y_val = [y[id] for id in idx_val]
	train_file = "training_data/dataset_features_ytrain.pkl.gz"
	val_file = "training_data/dataset_features_yval.pkl.gz"
	with gzip.open(train_file, "wb") as f:
		pickle.dump(y_train, f, protocol=4)
	with gzip.open(val_file, "wb") as f:
		pickle.dump(y_val, f, protocol=4)
	y = [train_file, val_file]


	m_train = [mol_order[id] for id in idx_train]
	m_val = [mol_order[id] for id in idx_val]
	train_file = "training_data/dataset_features_mtrain.pkl.gz"
	val_file = "training_data/dataset_features_mval.pkl.gz"
	with gzip.open(train_file, "wb") as f:
		pickle.dump(m_train, f, protocol=4)
	with gzip.open(val_file, "wb") as f:
		pickle.dump(m_val, f, protocol=4)
	mol_order = [train_file, val_file]

	return x, y, r, mol_order



# make BCAI features
def get_BCAI_features(mols, targetflag='CCS', training=True):

	target = flag_to_target(targetflag)

	atoms = BCAI.make_atom_df(mols)
	structure_dict = BCAI.make_structure_dict(atoms)
	BCAI.enhance_structure_dict(structure_dict)
	BCAI.enhance_atoms(atoms, structure_dict)

	bonds = BCAI.make_bonds_df(mols)
	bonds = BCAI.enhance_bonds(bonds, structure_dict)

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

	Dset = BCAI.create_dataset(atoms, bonds, triplets, labeled = True, max_count = 10**10, mol_order=mol_order)

	if training:
		x, y, r, mol_order = save_split_dataset(Dset, y, r, mol_order)
	else:
		x, y, r, mol_order = save_dataset(Dset, y, r, mol_order)

	return x, y, r, mol_order








































###

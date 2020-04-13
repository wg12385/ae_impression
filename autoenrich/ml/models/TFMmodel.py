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

import numpy as np
from autoenrich.ml.models.model import genericmodel
from sklearn.model_selection import KFold

import torch
torch.manual_seed(0)
from tqdm import tqdm
from torch.utils.data import TensorDataset, DataLoader
from torch import optim
from torch import autograd
import pickle
import gzip

import autoenrich.ml.models.BCAI.train as BCAI_train
import autoenrich.ml.models.BCAI.graph_transformer as BCAI_graph
import autoenrich.ml.models.BCAI.predictor as BCAI_predict
from .BCAI.modules import radam

import sys
import time
import random

class TFMmodel(genericmodel):

	def __init__(self, mol_order, id='TFMmodel', x=[], y=[], r=[], files=[], params={}, model_args={}):

		torch.cuda.empty_cache()

		genericmodel.__init__(self, id, x, y, params, model_args)
		self.r = r
		self.molfiles = files
		self.mol_order = mol_order

		self.N_atypes = [5, 10, 15]
		self.N_btypes = [30, 33, 33]
		self.N_ttypes = [23, 79]
		self.N_qtypes = [40]


		if len(x) > 0:
			with gzip.open(self.train_x[0], "rb") as f:
				train_dataset = TensorDataset(*pickle.load(f))
			self.N_atypes = [int(train_dataset.tensors[1][:,:,i].max()) for i in range(3)]   # Atom hierarchy has 3 levels
			self.N_btypes = [int(train_dataset.tensors[3][:,:,i].max()) for i in range(3)]   # Bond hierarchy has 3 levels
			self.N_ttypes = [int(train_dataset.tensors[5][:,:,i].max()) for i in range(2)]  # Triplet hierarchy has 2 levels
			self.N_qtypes = [int(train_dataset.tensors[7][:,:,i].max()) for i in range(1)]   # Quad hierarchy has only 1 level


	def train(self, train_x=[], train_y=[]):

		torch.backends.cudnn.deterministic = True
		torch.backends.cudnn.benchmark = False
		np.random.seed(0)
		random.seed(0)

		if len(train_x) == 0:
			train_x = self.train_x
		if len(train_y) == 0:
			train_y = self.train_y

		with gzip.open(train_x, "rb") as f:
			train_dataset = TensorDataset(*pickle.load(f))

		batch_size = 1

		if len(train_dataset[0]) <= batch_size:
			batch_size = len(train_dataset[0]) - 1

		train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=False)

		NUM_ATOM_TYPES = self.N_atypes  # Atom hierarchy has 3 levels
		NUM_BOND_TYPES = self.N_btypes   # Bond hierarchy has 3 levels
		NUM_TRIPLET_TYPES = self.N_ttypes  # Triplet hierarchy has 2 levels
		NUM_QUAD_TYPES = self.N_qtypes   # Quad hierarchy has only 1 level
		NUM_BOND_ORIG_TYPES = 8
		MAX_BOND_COUNT = 500  # params['max_bond_count']
		max_step = len(train_loader)

		device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

		d_model = int(self.params['d_model']/int(self.params['n_head'])*2)*int(self.params['n_head'])*2
		assert d_model % 2 == 0
		self.params['d_model'] = d_model
		self.model = BCAI_graph.GraphTransformer(dim=d_model, n_layers=int(self.params['n_layer']), d_inner=int(self.params['d_inner']),
								 fdim = 200, final_dim=int(self.params['final_dim']), dropout=self.params['dropout'],
								 dropatt=self.params['dropatt'], final_dropout=self.params['final_dropout'], n_head=int(self.params['n_head']),
								 num_atom_types=NUM_ATOM_TYPES,
								 num_bond_types=NUM_BOND_TYPES,
								 num_triplet_types=NUM_TRIPLET_TYPES,
								 dist_embedding='sine',
								 atom_angle_embedding='sine',
								 trip_angle_embedding='sine',
								 wnorm=True).to(device)

		if torch.cuda.device_count() > 1:
			print("Using", torch.cuda.device_count(), "GPUs")
			self.model = torch.nn.DataParallel(self.model)

		self.params['n_all_param'] = sum([p.nelement() for p in self.model.parameters() if p.requires_grad])

		self.params['optim'] = "RAdam"
		optimizer = getattr(optim if self.params['optim'] != "RAdam" else radam, self.params['optim'])(self.model.parameters(), lr=self.params['learning_rate'])
		scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, max_step, eta_min=self.params['eta_min'])

		para_model = self.model.to(device)

		pbar_iter = tqdm(range(int(self.params['tr_epochs'])), desc='epoch loss: ')
		for tr_epoch in pbar_iter:
			loss1, loss2, loss3, mabs = BCAI_train.epoch(train_loader, self.model, optimizer, self.params['learning_rate'])
			string = "epoch loss: {0:<10.4f} / {1:<10.4f}".format(mabs, loss1)
			pbar_iter.set_description(string)
		self.trained = True

	def predict(self, test_x):

		with gzip.open(test_x, "rb") as f:
			test_dataset = TensorDataset(*pickle.load(f))

		test_loader = DataLoader(test_dataset, batch_size=1, shuffle=True, drop_last=True)

		#BCAI_predict.single_model_predict(test_loader, self.model, 'test')
		MAX_BOND_COUNT = 500

		y_test, y_pred = BCAI_predict.single_model_predict(test_loader, self.model, "name")


		return y_test, y_pred

	def cv_predict(self, fold):

		#with autograd.detect_anomaly():
		if fold != 1:
			print('For TFM models, single fold CV validation is always performed. sorry.')

		molnames = self.mol_order

		self.train(train_x=self.train_x[0], train_y=self.train_y[0])

		test, preds = self.predict(self.train_x[1])

		pred_y = np.asarray(preds)

		return np.mean(np.absolute(pred_y - np.asarray(test)))



















##

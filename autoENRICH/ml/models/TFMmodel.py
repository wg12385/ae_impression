# Copyright 2020 Will Gerrard
#This file is part of autoENRICH.

#autoENRICH is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoENRICH is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoENRICH.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from autoENRICH.ml.models.model import genericmodel
from sklearn.model_selection import KFold

import torch
from tqdm import tqdm
from torch.utils.data import TensorDataset, DataLoader
from torch import optim
from torch import autograd
import pickle
import gzip

import autoENRICH.ml.models.BCAI.train as BCAI_train
import autoENRICH.ml.models.BCAI.graph_transformer as BCAI_graph
import autoENRICH.ml.models.BCAI.predictor as BCAI_predict
from .BCAI.modules import radam

import sys

import random

class TFMmodel(genericmodel):


	def __init__(self, mol_order, id='TFMmodel', x=[], y=[], r=[], files=[], params={}, model_args={}):
		genericmodel.__init__(self, id, x, y, params, model_args)
		self.r = r
		self.molfiles = files
		self.mol_order = mol_order

		self.N_atypes = [5, 10, 15]
		self.N_btypes = [30, 33, 33]
		self.N_ttypes = [23, 79]
		self.N_qtypes = [40]


		if len(x) > 0:
			train_dataset = TensorDataset(*self.train_x)
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

		train_dataset = TensorDataset(*train_x)
		train_loader = DataLoader(train_dataset, batch_size=10, shuffle=True, drop_last=True)

		NUM_ATOM_TYPES = self.N_atypes  # Atom hierarchy has 3 levels
		NUM_BOND_TYPES = self.N_btypes   # Bond hierarchy has 3 levels
		NUM_TRIPLET_TYPES = self.N_ttypes  # Triplet hierarchy has 2 levels
		NUM_QUAD_TYPES = self.N_qtypes   # Quad hierarchy has only 1 level
		NUM_BOND_ORIG_TYPES = 8
		MAX_BOND_COUNT = 500  # params['max_bond_count']
		max_step = len(train_loader)

		device = torch.device('cuda')

		self.model = BCAI_graph.GraphTransformer(dim=200, n_layers=int(self.params['n_layer']), d_inner=600,
								 fdim = 200, final_dim=int(self.params['final_dim']), dropout=self.params['dropout'],
								 dropatt=self.params['dropatt'], final_dropout=self.params['final_dropout'], n_head=10,
								 num_atom_types=NUM_ATOM_TYPES,
								 num_bond_types=NUM_BOND_TYPES,
								 num_triplet_types=NUM_TRIPLET_TYPES,
								 num_quad_types=NUM_QUAD_TYPES,
								 dist_embedding='sine',
								 atom_angle_embedding='sine',
								 trip_angle_embedding='sine',
								 quad_angle_embedding='sine',
								 wnorm=True,
								 use_quad=False).to(device)

		self.params['n_all_param'] = sum([p.nelement() for p in self.model.parameters() if p.requires_grad])

		self.params['optim'] = "RAdam"
		optimizer = getattr(optim if self.params['optim'] != "RAdam" else radam, self.params['optim'])(self.model.parameters(), lr=self.params['learning_rate'])
		scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, max_step, eta_min=self.params['eta_min'])

		para_model = self.model.to(device)
		train_epochs = 10
		for tr_epoch in range(train_epochs):
			print('\ttrepoch: ', tr_epoch, '/', train_epochs)
			loss1, loss2, loss2 = BCAI_train.epoch(train_loader, self.model, optimizer, self.params['learning_rate'])
		self.trained = True
		'''
		print('PARAMS:')
		with open('named_model_params.txt', 'w') as f:
			for name, param in self.model.named_parameters():
				print('PARAMETER ---------------------------', file=f)
				print(name, file=f)
				print(param.data, file=f)

		print(self.model)
		'''

	def predict(self, test_x, train_x=[]):

		test_dataset = TensorDataset(*test_x)
		test_loader = DataLoader(test_dataset, batch_size=1, shuffle=True, drop_last=True)

		#BCAI_predict.single_model_predict(test_loader, self.model, 'test')
		MAX_BOND_COUNT = 500
		dev = "cuda"
		self.model.to(dev)
		self.model.eval()
		y_predictions = []
		with torch.no_grad():
			for arr in tqdm(test_loader):
				x_idx, x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, x_quad, x_quad_angle, y = arr
				x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, x_quad, x_quad_angle, y = \
					x_atom.to(dev), x_atom_pos.to(dev), x_bond.to(dev), x_bond_dist.to(dev), \
					x_triplet.to(dev), x_triplet_angle.to(dev), x_quad.to(dev), x_quad_angle.to(dev), y.to(dev)

				x_bond, x_bond_dist, y = x_bond[:, :MAX_BOND_COUNT], x_bond_dist[:, :MAX_BOND_COUNT], y[:,:MAX_BOND_COUNT]
				y_pred, _ = self.model(x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, x_quad, x_quad_angle)
				y_pred_pad = torch.cat([torch.zeros(y_pred.shape[0], 1, y_pred.shape[2], device=y_pred.device), y_pred], dim=1)
				y_pred_scaled = y_pred_pad.gather(1,x_bond[:,:,1][:,None,:])[:,0,:] * y[:,:,2] + y[:,:,1]

				y_selected = y_pred_scaled.masked_select((x_bond[:,:,0] > 0) & (y[:,:,3] > 0)).cpu().numpy()
				ids_selected = y[:,:,0].masked_select((x_bond[:,:,0] > 0) & (y[:,:,3] > 0))
				ids_selected = ids_selected.numpy()

				for id_, pred in zip(ids_selected, y_selected):
					y_predictions.append(pred)

		return y_predictions

	def cv_predict(self, fold):

		#with autograd.detect_anomaly():

			if fold >= 2:
				molnames = self.mol_order

				kf = KFold(n_splits=fold)
				kf.get_n_splits(self.train_x[0])
				pred_y = []

				train_dataset = TensorDataset(*self.train_x)
				self.N_atypes = [int(train_dataset.tensors[1][:,:,i].max()) for i in range(3)]   # Atom hierarchy has 3 levels
				self.N_btypes = [int(train_dataset.tensors[3][:,:,i].max()) for i in range(3)]   # Bond hierarchy has 3 levels
				self.N_ttypes = [int(train_dataset.tensors[5][:,:,i].max()) for i in range(2)]  # Triplet hierarchy has 2 levels
				self.N_qtypes = [int(train_dataset.tensors[7][:,:,i].max()) for i in range(1)]   # Quad hierarchy has only 1 level


				for train_index, test_index in kf.split(self.train_x[0]):

					train_x_list = []
					train_y_list = []
					test_x_list = []
					test_y_list = []

					for i in range(len(self.train_x)):
						train_x_list.append(torch.index_select(self.train_x[i], 0, torch.tensor(train_index)))
						test_x_list.append(torch.index_select(self.train_x[i], 0, torch.tensor(test_index)))

					for r, ref in enumerate(self.r):
						if ref[0] in [molnames[idx] for idx in train_index]:
							train_y_list.append(self.train_y[r])
						else:
							test_y_list.append(self.train_y[r])

					assert len(train_x_list) == len(train_x_list)

					self.train(train_x=train_x_list, train_y=train_y_list)
					preds = self.predict(test_x_list)
					pred_y.extend(self.predict(test_x_list))

			elif fold == 1:
				# Splits as in 5-fold, but just does 1 permutation
				molnames = self.mol_order

				kf = KFold(n_splits=5)
				kf.get_n_splits(self.train_x[0])
				pred_y = []

				train_dataset = TensorDataset(*self.train_x)
				self.N_atypes = [int(train_dataset.tensors[1][:,:,i].max()) for i in range(3)]   # Atom hierarchy has 3 levels
				self.N_btypes = [int(train_dataset.tensors[3][:,:,i].max()) for i in range(3)]   # Bond hierarchy has 3 levels
				self.N_ttypes = [int(train_dataset.tensors[5][:,:,i].max()) for i in range(2)]  # Triplet hierarchy has 2 levels
				self.N_qtypes = [int(train_dataset.tensors[7][:,:,i].max()) for i in range(1)]   # Quad hierarchy has only 1 level

				# This is a bad way to do this, needs fixing
				# Priority to make sure same as multiple fold atm
				for train_index, test_index in kf.split(self.train_x[0]):

					train_x_list = []
					train_y_list = []
					test_x_list = []
					test_y_list = []


					for i in range(len(self.train_x)):
						train_x_list.append(torch.index_select(self.train_x[i], 0, torch.tensor(train_index)))
						test_x_list.append(torch.index_select(self.train_x[i], 0, torch.tensor(test_index)))

					for r, ref in enumerate(self.r):
						if ref[0] in [molnames[idx] for idx in train_index]:
							train_y_list.append(self.train_y[r])
						else:
							test_y_list.append(self.train_y[r])

					self.train(train_x=train_x_list, train_y=train_y_list)
					#preds = self.predict(test_x_list)
					pred_y.extend(self.predict(test_x_list))

					pred_y = np.asarray(pred_y)

					print(pred_y)
					print(test_y_list)

					return np.mean(np.absolute(pred_y - np.asarray(test_y_list)))

			else:
				train_dataset = TensorDataset(*self.train_x)
				self.N_atypes = [int(train_dataset.tensors[1][:,:,i].max()) for i in range(3)]   # Atom hierarchy has 3 levels
				self.N_btypes = [int(train_dataset.tensors[3][:,:,i].max()) for i in range(3)]   # Bond hierarchy has 3 levels
				self.N_ttypes = [int(train_dataset.tensors[5][:,:,i].max()) for i in range(2)]  # Triplet hierarchy has 2 levels
				self.N_qtypes = [int(train_dataset.tensors[7][:,:,i].max()) for i in range(1)]   # Quad hierarchy has only 1 level

				pred_y = []
				self.train(train_x=self.train_x, train_y=self.train_y)
				pred_y.extend(self.train_y)

			pred_y = np.asarray(pred_y)

			return pred_y



















##

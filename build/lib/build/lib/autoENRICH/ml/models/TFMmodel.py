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
import pickle
import gzip

import autoENRICH.ml.models.BCAI.train as BCAI_train
import autoENRICH.ml.models.BCAI.graph_transformer as BCAI_graph
import autoENRICH.ml.models.BCAI.predictor as BCAI_predict
from .BCAI.modules import radam

import sys

class TFMmodel(genericmodel):


	def __init__(self, id='TFMmodel', x=[], y=[], r=[], files=[], params={}, model_args={}):
		genericmodel.__init__(self, id, x, y, params, model_args)
		self.r = r
		self.molfiles = files


	def train(self, train_x=[], train_y=[]):

		if len(train_x) == 0:
			train_x = self.train_x
		if len(train_y) == 0:
			train_y = self.train_y

		train_dataset = TensorDataset(*train_x)
		train_loader = DataLoader(train_dataset, batch_size=10, shuffle=True, drop_last=True)

		NUM_ATOM_TYPES = [int(train_dataset.tensors[1][:,:,i].max()) for i in range(3)]   # Atom hierarchy has 3 levels
		NUM_BOND_TYPES = [int(train_dataset.tensors[3][:,:,i].max()) for i in range(3)]   # Bond hierarchy has 3 levels
		NUM_TRIPLET_TYPES = [int(train_dataset.tensors[5][:,:,i].max()) for i in range(2)]  # Triplet hierarchy has 2 levels
		NUM_QUAD_TYPES = [int(train_dataset.tensors[7][:,:,i].max()) for i in range(1)]   # Quad hierarchy has only 1 level
		NUM_BOND_ORIG_TYPES = 8
		MAX_BOND_COUNT = 500  # params['max_bond_count']
		max_step = len(train_loader)

		device = torch.device('cpu')

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

		self.params['optim'] = "Adam"
		optimizer = getattr(optim if self.params['optim'] != "RAdam" else radam, self.params['optim'])(self.model.parameters(), lr=self.params['learning_rate'])
		scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, max_step, eta_min=self.params['eta_min'])

		para_model = self.model.to(device)

		BCAI_train.epoch(train_loader, self.model, optimizer, self.params['learning_rate'])
		self.trained = True


	def predict(self, test_x, train_x=[]):

		test_dataset = TensorDataset(*test_x)
		test_loader = DataLoader(test_dataset, batch_size=1, shuffle=True, drop_last=True)

		#BCAI_predict.single_model_predict(test_loader, self.model, 'test')
		MAX_BOND_COUNT = 500
		dev = "cpu"
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

		molnames = list(set([self.r[i][0] for i in range(len(self.r))]))

		kf = KFold(n_splits=fold)
		kf.get_n_splits(self.train_x[0])
		pred_y = []

		for train_index, test_index in kf.split(self.train_x[0]):

			train_x_list = []
			train_y_list = []
			test_x_list = []
			test_y_list = []


			list1 = [self.train_x[0][idx] for idx in train_index]
			tensor1 = torch.tensor(list1)


			for i in range(len(self.train_x)):
				train_x_list.append(torch.index_select(self.train_x[i], 0, torch.tensor(train_index)))
				test_x_list.append(torch.index_select(self.train_x[i], 0, torch.tensor(test_index)))

			for r, ref in enumerate(self.r):
				if ref[0] in [molnames[idx] for idx in train_index]:
					train_y_list.append(self.train_y[r])
				else:
					test_y_list.append(self.train_y[r])

			self.train(train_x=train_x_list, train_y=train_y_list)
			pred_y.extend(self.predict(test_x_list))

		pred_y = np.asarray(pred_y)

		return pred_y



















##
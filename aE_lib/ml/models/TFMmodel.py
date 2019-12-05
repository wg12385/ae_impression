import numpy as np
from ml.models.model import genericmodel
from sklearn.model_selection import KFold

import torch
from tqdm import tqdm
from torch.utils.data import TensorDataset, DataLoader


import .BCAI.train as BCAI_train
import .BCAI.graph_transformer as BCAI_graph
import .BCAI.predictor as BCAI_predict

class TFMmodel(genericmodel):


	def __init__(self):
		genericmodel.__init__(self, id, x, y, params, model_args)



	def train(self):

		train_dataset = TensorDataset(x)
		train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=True)

		NUM_ATOM_TYPES = [int(train_dataset.tensors[1][:,:,i].max()) for i in range(3)]   # Atom hierarchy has 3 levels
		NUM_BOND_TYPES = [int(train_dataset.tensors[3][:,:,i].max()) for i in range(3)]   # Bond hierarchy has 3 levels
		NUM_TRIPLET_TYPES = [int(train_dataset.tensors[5][:,:,i].max()) for i in range(2)]  # Triplet hierarchy has 2 levels
		NUM_QUAD_TYPES = [int(train_dataset.tensors[7][:,:,i].max()) for i in range(1)]   # Quad hierarchy has only 1 level
		NUM_BOND_ORIG_TYPES = 8
		MAX_BOND_COUNT = 500  # params['max_bond_count']
		max_step = len(train_loader)

		device = torch.device('cpu')
		self.model = BCAI_graph.GraphTransformer(dim=int(params['d_model']), n_layers=int(params['n_layer']), d_inner=3800,
								 fdim = 200, final_dim=int(params['final_dim']), dropout=params['dropout'],
								 dropatt=params['dropatt'], final_dropout=params['final_dropout'], n_head=int(params['n_head']),
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
		params['n_all_param'] = sum([p.nelement() for p in model.parameters() if p.requires_grad])

		params['optim'] = "RAdam"
		optimizer = getattr(optim if params['optim'] != "RAdam" else radam, params['optim'])(self.model.parameters(), lr=params['lr'])
		scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, max_step, eta_min=params['eta_min'])

		para_model = self.model.to(device)

		BCAI_train.epoch(train_loader, self.model, optimizer)
		self.trained = True


	def predict(self, test_x):

		test_dataset = TensorDataset(test_x)
		test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=True, drop_last=True)

		#BCAI_predict.single_model_predict(test_loader, self.model, 'test')
		MAX_BOND_COUNT = 500
		#dev = "cpu"
		model = model.to(dev)
		model.eval()
		with torch.no_grad():
			for arr in tqdm(loader):
				x_idx, x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, x_quad, x_quad_angle, y = arr
				x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, x_quad, x_quad_angle, y = \
					x_atom.to(dev), x_atom_pos.to(dev), x_bond.to(dev), x_bond_dist.to(dev), \
					x_triplet.to(dev), x_triplet_angle.to(dev), x_quad.to(dev), x_quad_angle.to(dev), y.to(dev)

				x_bond, x_bond_dist, y = x_bond[:, :MAX_BOND_COUNT], x_bond_dist[:, :MAX_BOND_COUNT], y[:,:MAX_BOND_COUNT]
				y_pred, _ = model(x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, x_quad, x_quad_angle)

		return y_pred




















##

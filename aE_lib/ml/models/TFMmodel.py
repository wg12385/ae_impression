import numpy as np
from ml.models.model import genericmodel
from sklearn.model_selection import KFold

import torch
from tqdm import tqdm
from torch.utils.data import TensorDataset, DataLoader
from torch import optim
import pickle
import gzip

import ml.models.BCAI.train as BCAI_train
import ml.models.BCAI.graph_transformer as BCAI_graph
import ml.models.BCAI.predictor as BCAI_predict
from .BCAI.modules import radam

class TFMmodel(genericmodel):


	def __init__(self, id='TFMmodel', x=[], y=[], params={}, model_args={}):
		genericmodel.__init__(self, id, x, y, params, model_args)



	def train(self, train_x=[], train_y=[]):

		if len(train_x) == 0:
			train_x = self.train_x
		if len(train_y) == 0:
			train_y = self.train_y

		with gzip.open("torch_proc_submission.pkl.gz", "rb") as f:
			train_dataset = TensorDataset(*pickle.load(f))
		train_loader = DataLoader(train_dataset, batch_size=10, shuffle=True, drop_last=True)

		NUM_ATOM_TYPES = [int(train_dataset.tensors[1][:,:,i].max()) for i in range(3)]   # Atom hierarchy has 3 levels
		NUM_BOND_TYPES = [int(train_dataset.tensors[3][:,:,i].max()) for i in range(3)]   # Bond hierarchy has 3 levels
		NUM_TRIPLET_TYPES = [int(train_dataset.tensors[5][:,:,i].max()) for i in range(2)]  # Triplet hierarchy has 2 levels
		NUM_QUAD_TYPES = [int(train_dataset.tensors[7][:,:,i].max()) for i in range(1)]   # Quad hierarchy has only 1 level
		NUM_BOND_ORIG_TYPES = 8
		MAX_BOND_COUNT = 500  # params['max_bond_count']
		max_step = len(train_loader)

		device = torch.device('cpu')


		self.model = BCAI_graph.GraphTransformer(dim=200, n_layers=int(self.params['n_layer']), d_inner=3800,
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

		'''
		self.model = BCAI_graph.GraphTransformer(dim=int(self.params['d_model']), n_layers=int(self.params['n_layer']), d_inner=3800,
								 fdim = 200, final_dim=int(self.params['final_dim']), dropout=self.params['dropout'],
								 dropatt=self.params['dropatt'], final_dropout=self.params['final_dropout'], n_head=int(self.params['n_head']),
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
		'''
		self.params['n_all_param'] = sum([p.nelement() for p in self.model.parameters() if p.requires_grad])

		self.params['optim'] = "Adam"
		optimizer = getattr(optim if self.params['optim'] != "RAdam" else radam, self.params['optim'])(self.model.parameters(), lr=self.params['learning_rate'])
		scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, max_step, eta_min=self.params['eta_min'])

		para_model = self.model.to(device)

		BCAI_train.epoch(train_loader, self.model, optimizer, self.params['learning_rate'])
		self.trained = True


	def predict(self, test_x, train_x=[]):

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

	def cv_predict(self, fold):



		'''
		kf = KFold(n_splits=fold)

		kf.get_n_splits(len(self.train_x[0]))
		pred_y = []

		for train_index, test_index in kf.split(self.train_x):
			train_x = []
			test_x = []
			for tr_x in self.train_x:
				train_x.append(torch.as_tensor(np.asarray(tr_x)[train_index]))
				test_x.append(torch.as_tensor(np.asarray(tr_x)[test_index]))


			self.train(train_x=train_x, train_y=[])

			pred_y.extend(self.predict(test_x, train_x=[]))
		'''

		self.train(train_x=self.train_x)
		pred_y = self.predict(self.train_x)

		pred_y = np.asarray(pred_y)

		return pred_y


















##

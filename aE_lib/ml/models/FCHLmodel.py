import numpy as np
from ml.models.model import genericmodel
from sklearn.model_selection import KFold
import qml
import copy

class FCHLmodel(genericmodel):

	def __init__(self, id='FCHLmodel', x=[], y=[], params={}, model_args={}):
		genericmodel.__init__(self, id, x, y, params, model_args)


	def train(self):
		k = []

		dimensions = self.train_x.shape[1]

		# reshape x array:
		Xr = []
		for _ in range(dimensions):
			Xr.append([])
		for x in self.train_x:
			for i in range(dimensions):
				Xr[i].append(x[i])
		Xr = np.asarray(Xr)

		# loop through representation sets in training set
		for d in range(dimensions):
			# append kernel for each representation set
			k.append(qml.fchl.get_atomic_symmetric_kernels(Xr[d], verbose=False,
					two_body_scaling=self.params['two_body_scaling'], three_body_scaling=self.params['three_body_scaling'],
					two_body_width=self.params['two_body_width'], three_body_width=self.params['three_body_width'],
					two_body_power=self.params['two_body_power'], three_body_power=self.params['three_body_power'],
					cut_start=self.params['cut_start'], cut_distance=self.params['cutoff'],
					fourier_order=1, alchemy="periodic-table",
					alchemy_period_width=self.params['alchemy_period_width'], alchemy_group_width=self.params['alchemy_group_width'],
					kernel="gaussian", kernel_args={'sigma':[self.params['sigma']]})[0] + self.params['lamda'] * np.identity(Xr[d].shape[0]))
		# loop through kernels
		K = k[0]
		if len(k) > 1:
			for i in range(1, len(k)):
				# multiply kernels so result is k1 * k2 * k3 ...
				K = K * k[i]

		# get the KRR prefactors, i.e. this is the training of the network
		self.alpha = qml.math.cho_solve(K, self.train_y)

		# report training state
		self.trained = True


	def predict(self, test_x):
		ks = []
		assert test_x.shape[1] == self.train_x.shape[1]

		dimensions = self.train_x.shape[1]

		# reshape x arrays:
		Xe = []
		Xr = []
		for _ in range(dimensions):
			Xe.append([])
			Xr.append([])
		for x in test_x:
			for i in range(dimensions):
				Xe[i].append(x[i])
		for x in self.train_x:
			for i in range(dimensions):
				Xr[i].append(x[i])
		Xe = np.asarray(Xe)
		Xr = np.asarray(Xr)

		# loop through representation sets in training set
		for d in range(dimensions):
			# append kernel for each representation set
			ks.append(qml.fchl.get_atomic_kernels(Xe[d], Xr[d], verbose=False, \
					two_body_scaling=np.sqrt(8), three_body_scaling=1.6,
					two_body_width=0.2, three_body_width=np.pi,
					two_body_power=4.0, three_body_power=2.0,
					cut_start=1.0, cut_distance=self.params['cutoff'],
					fourier_order=1, alchemy="periodic-table",
					alchemy_period_width=1.6, alchemy_group_width=1.6,
					kernel="gaussian", kernel_args={'sigma':[self.params['sigma']]})[0])


		# loop through kernels
		Ks = ks[0]
		if len(ks) > 1:
			for i in range(1, len(ks)):
				# multiply kernels so result is k1 * k2 * k3 ...
				Ks = Ks * ks[i]

		# predict values of y
		y_pred = np.dot(Ks, self.alpha)

		return y_pred

	def cv_predict(self, fold):
		kf = KFold(n_splits=fold)
		kf.get_n_splits(self.train_x)
		pred_y = []
		for train_index, test_index in kf.split(self.train_x):

			tmp_train_x = np.asarray(self.train_x[train_index])
			tmp_train_y = self.train_y[train_index]
			tmp_test_x = np.asarray(self.train_x[test_index])
			tmp_test_y = self.train_y[test_index]

			tmp_model = copy.deepcopy(self)
			tmp_model.train_x = tmp_train_x
			tmp_model.train_y = tmp_train_y

			tmp_model.train()

			pred_y.extend(tmp_model.predict(tmp_test_x))

		pred_y = np.asarray(pred_y)

		return pred_y

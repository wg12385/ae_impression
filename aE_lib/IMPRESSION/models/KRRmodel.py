
from IMPRESSION.models.model import genericmodel
import qml

class KRRmodel(genericmodel):

	def __init__(self, id='KRRmodel', x=[], y=[], params={}):
		genericmodel.__init__(self, id, x, y, params)
		self.alpha = []


	def train(self):
		k = []
		# loop through representation sets in training set

		# reshape x array:
		Xr = []
		for _ in dimensions:
			Xr.append([])
		for x in self.train_x:
			for i in range(dimensions):
				Xr[i].append(x[i])

		# loop through representation sets in training set
		for d in range(dimensions):
			# append kernel for each representation set
			k.append(qml.kernels.laplacian_kernel(Xr[d], Xr[d], self.params['sig']) + self.params['lam'] * np.identity(Xr[d].shape[0]))

		# loop through kernels
		K = k[0]
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

		dimensions = train_x.shape[1]

		# reshape x arrays:
		Xe = []
		Xr = []
		for _ in dimensions:
			Xe.append([])
			Xr.append([])
		for x in self.test_x:
			for i in range(dimensions):
				Xe[i].append(x[i])
		for x in self.train_x:
			for i in range(dimensions):
				Xr[i].append(x[i])

		# loop through representation sets in training set
		for d in range(dimensions):
			# append kernel for each representation set
			ks.append(qml.kernels.laplacian_kernel(Xe[d], Xr[d], self.params['sig']))

		# loop through kernels
		Ks = ks[0]
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

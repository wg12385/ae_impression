import qml
import numpy as np
import copy
from qml.math import cho_solve
from qml.kernels import laplacian_kernel
from sklearn.model_selection import KFold
# top level model class

class genericmodel(object):

	def __init__(self, id='genericmodel', train_x=[], train_y=[], params={})

		self.train_x = train_x
		self.train_y = train_y

		self.params = params



	def cv_predict(self, fold):
		kf = KFold(n_splits=fold)
		kf.get_n_splits(self.train_x[0])
		pred_y = []
		for train_index, test_index in kf.split(self.train_x[0]):
			sets = []
			train_sets = []
			test_sets = []
			for set in self.train_x:
				train_sets.append(set[train_index])
				test_sets.append(set[test_index])

			tmp_train_x = np.asarray(train_sets)
			tmp_train_y = self.train_y[train_index]
			tmp_test_x = np.asarray(test_sets)
			tmp_test_y = self.train_y[test_index]

			tmp_model = copy.deepcopy(self)
			tmp_model.train_x = tmp_train_x
			tmp_model.train_y = tmp_train_y

			tmp_model.train()

			pred_y.extend(tmp_model.predict(tmp_test_x))

		pred_y = np.asarray(pred_y)

		return pred_y

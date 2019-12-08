
import qml

import numpy as np
import copy
from qml.math import cho_solve
from qml.kernels import laplacian_kernel
from sklearn.model_selection import KFold

import tracemalloc
# top level model class

class genericmodel(object):

	def __init__(self, id='genericmodel', train_x=[], train_y=[], params={},
									model_args={}):

		self.train_x = train_x
		self.train_y = train_y

		self.params = params
		self.args = model_args

		self.trained = False


	def train(self):
		print('Stub function, why are you running this ??')

	def predict(self):
		print('Stub function, why are you running this ??')


	def cv_predict(self, fold):
		kf = KFold(n_splits=fold)
		kf.get_n_splits(np.asarray(self.train_x))
		pred_y = []

		for train_index, test_index in kf.split(np.asarray(self.train_x)):

			self.train(train_x=np.asarray(self.train_x)[train_index], train_y=np.asarray(self.train_y)[train_index])

			pred_y.extend(self.predict(np.asarray(self.train_x)[test_index], train_x=np.asarray(self.train_x)[train_index]))

		pred_y = np.asarray(pred_y)

		return pred_y


























##

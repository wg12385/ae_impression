import numpy as np
import copy

from ml.models.model import genericmodel
from sklearn.model_selection import KFold
import sys

import tensorflow as tf
from keras.models import Sequential
from keras.layers import Dense, Activation, Flatten
from keras.optimizers import RMSprop
from keras.callbacks import EarlyStopping


class NNmodel(genericmodel):
	def __init__(self, id='NNmodel', x=[], y=[], params={}, model_args={}):
		genericmodel.__init__(self, id, x, y, params, model_args)



	def train(self):
		self.net = Sequential()

		print(self.train_x.shape)
		print('FIX CODE FOR THIS SHAPE')

		dimensions = self.train_x.shape[1]

		# reshape x array:
		Xr = []
		for _ in range(dimensions):
			Xr.append([])
		for x in self.train_x:
			for i in range(dimensions):
				Xr[i].extend(x[i])
		Xr = np.asarray(Xr)


		self.net.add(Dense(self.params['hidden_neurons'], input_shape=Xr.shape, kernel_initializer='random_uniform'))
		self.net.add(Activation('relu'))

		if self.params['hidden_layers'] > 1:
			for layer in range(self.params['hidden_layers']-1):
				self.net.add(Dense(self.params['hidden_neurons'], input_shape=Xr.shape, kernel_initializer='random_uniform'))
				self.net.add(Activation('relu'))
		self.net.add(Flatten())

		self.net.add(Dense(1, input_dim=self.params['hidden_neurons'], kernel_initializer='random_uniform'))
		optzer = RMSprop(lr=10**self.params['learning_rate'], rho=0.9, epsilon=None, decay=0.0)
		self.net.compile(optimizer=optzer, loss='mse')

		es = EarlyStopping(monitor='loss', mode='min', verbose=1)

		self.net.fit(Xr, self.train_y, epochs=self.params['nn_epochs'], batch_size=self.params['batch_size'], verbose=1, callbacks=[es])
		self.trained = True

	def predict(self, test_x):

		dimensions = test_x.shape[1]

		# reshape x array:
		Xe = []
		for _ in range(dimensions):
			Xe.append([])
		for x in test_x:
			for i in range(dimensions):
				Xe[i].extend(x[i])
		Xe = np.asarray(Xe)


		pred_y = self.net.predict(Xe, batch_size=None)


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















###

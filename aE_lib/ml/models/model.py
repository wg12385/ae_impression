
import qml

import numpy as np
import copy
from qml.math import cho_solve
from qml.kernels import laplacian_kernel
from sklearn.model_selection import KFold
# top level model class

class genericmodel(object):

	def __init__(self, id='genericmodel', train_x=[], train_y=[], params={},
									model_args={}):

		self.train_x = train_x
		self.train_y = train_y

		self.params = params
		self.args = model_args

		self.trained = False

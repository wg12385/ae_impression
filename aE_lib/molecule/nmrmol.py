import numpy as np
from file_read import orca_read, g09_read, structure_read
import pickle


class nmrmol(object):
	"""
		nmr data file
	"""

	def __init__(self, molid, path=''):
		self.path = path

		self.molid = str(molid)
		self.types = []
		self.xyz = []
		self.conn = []

		self.shift = []
		self.shift_var = []

		self.coupling = []
		self.coupling_var = []
		self.coupling_len = []

		self.energy = -404.404


	def read_structure(self, file, type):
		old_type_array = self.types
		old_xyz_array = self.xyz
		self.xyz, self.types, self.conn, self.coupling_len = structure_read.generic_pybel_read(file, type)

		if not (np.array_equal(old_type_array, self.types) or np.array_equal(old_xyz_array, self.xyz)):
			self.shift = []
			self.shift_var = []
			self.coupling = []
			self.coupling_var = []

	def read_opt(self, file, type):
		if type == 'orca':
			#self.xyz, self.types, self.conn, self.coupling_len = orca_read.read_structure(file)
			self.energy = orca_read.read_opt(file)


	def read_nmr(self, file, type):
		if type == 'orca':
			self.xyz, self.types, self.conn, self.coupling_len = orca_read.read_structure(file)
			self.shift, self.coupling = orca_read.read_nmr(file)
		if type == 'g09':
			self.xyz, self.types, self.conn, self.coupling_len = structure_read.generic_pybel_read(file, type)
			self.shift, self.coupling = g09_read.read_nmr(file)
		if type == 'nmredata':
			self.xyz, self.types, self.conn, self.coupling_len = structure_read.generic_pybel_read(file, 'sdf')
			self.shift, self.shift_var, self.coupling, self.coupling_var = nmredata_read.read_nmr(file)

	def save_pickle(file):
		pickle.dump(self, open(file, "wb"))














####

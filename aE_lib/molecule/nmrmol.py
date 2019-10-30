import numpy as np
from file_read import orca_read, g09_read, nmredata_read, structure_read
import pickle


class nmrmol(object):
	"""
		nmr data file
	"""

	def __init__(self, name, path='', molid=0):
		self.name = name
		self.path = path

		self.molid = molid
		self.types = np.zeros(len(5), dtype=np.int32)
		self.xyz = np.zeros((len(self.types),3) dtype=np.float64)
		self.conn = np.zeros((len(self.types), len(self.types)), dtype=np.float32)

		self.shift = np.zeros(len(self.types), dtype=np.float64)
		self.shift_var = np.zeros(len(self.types), dtype=np.float64)

		self.coupling = np.zeros((len(self.types), len(self.types)), dtype=np.float64)
		self.coupling_var = np.zeros((len(self.types), len(self.types)), dtype=np.float64)
		self.coupling_len = np.zeros((len(self.types), len(self.types)), dtype=np.int32)

		self.energy = -404.404


	def read_structure(file, type):
		old_type_array = self.types
		old_xyz_array = self.xyz
		self.xyz, self.types, self.conn, self.coupling_len = structure_read.generic_pybel_read(file, type)

		if not ( np.array_equal(old_type_array, self.types) or np.array_equal(old_xyz_array, self.xyz)):
			self.shift = np.zeros(len(types), dtype=np.float64)
			self.shift_var = np.zeros(len(types), dtype=np.float64)
			self.coupling = np.zeros((len(types), len(types)), dtype=np.float64)
			self.coupling_var = np.zeros((len(types), len(types)), dtype=np.float64)

	def read_opt(file, type):
		if type == 'orca':
			self.xyz, self.types, self.conn, self.coupling_len = orca_read.read_structure(file)
			self.energy = orca_read.read_opt(file)


	def read_nmr(file, type):
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

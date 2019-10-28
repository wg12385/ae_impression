import numpy as np
from conversion.convert import *

class nmrmol(object):
	"""
		nmr data file
	"""


	def __init__(self, name, path='', type='nxyz', pybmol='none', molid='', dummy='True'):
		self.name = name
		self.path = path

		if not dummy:
			if type == 'nxyz':
				self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = nxyz_to_nmrdata(name, path)

			elif type == 'g09':
				self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = g09_to_nmrdata(str(name+path))

			elif type == 'MIN_g09':
				self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = MINIMIZE_g09_to_nmrdata(str(name+path))

			elif type == 'tensor_g09':
				self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = TENSOR_g09_to_nmrdata(str(name+path))

			elif type == 'sygcml':
				self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = sygcml_to_nmrdata(name, path)

			elif type == 'pdb':
				self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = pdb_to_nmrdummy(str(name+path))

			elif type == 'mol2':
				self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = mol2_to_nmrdummy(str(name+path))

			elif type == 'xyz':
				self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = xyz_to_nmrdummy(str(name+path))

			elif type == 'pybmol':
				if pybmol == 'none':
					print('Must provid pybel mol object')
					sys.exit(0)
				else:
					self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = pybmol_to_nmrmol(pybmol, molid)

			else:
				print('TYPE NOT RECOGNISED:', type )
				sys.exit(0)

		if dummy:
			self.molid, self.types, self.xyz, self.dist = generic_pybel_read(file, type)
			self.shift = np.zeros(len(types), dtype=np.float64)
			self.shift_var = np.zeros(len(types), dtype=np.float64)
			self.coupling1b = []
			self.var1b = []
			self.coupling2b = []
			self.var2b = []
			self.coupling3b = []
			self.var3b = []


		self.molid = str(self.molid).lower()

	def make_pybmol(self):
		file = self.molid + 'tmp.xyz'

		nmrmol_to_xyz(self, file)

		pybmol = next(pyb.readfile('xyz', file))

		os.remove(file)

		return pybmol, self.molid







####

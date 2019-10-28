import numpy as np
import glob
from conversion.convert import *
from conversion.read_orca import eread, check_opt_status, check_nmr_status
from conversion.make_orca import make_optin, make_nmrin, make_submission_array, make_submission_qsub

class conformer(object):

	def __init__(self, xyz, types, id):

		# essential parameters
		self.xyz = xyz
		self.types = types
		self.molid = id

		# placeholders, optional at initialisation
		self.shift = np.zeros(len(types), dtype=np.float64)
		self.coupling1b = []
		self.coupling2b = []
		self.coupling3b = []
		self.coupling4b = []
		self.coupling5b = []

		# store location of conf search xyz file
		self.xyz_file = 'None'

		# store location and status of optimisation files
		self.opt_in = 'None'
		self.opt_log = 'None'
		self.opt_status = 'None'

		# store location and status of NMR files
		self.nmr_in = 'None'
		self.nmr_log = 'None'
		self.nmr_status = 'None'

		self.energy = 404.40404
		self.pop = 4.0404

	def get_energy_from_file(self, file):
		self.energy = eread(file)

	def data_from_file(self, file, type='g09', molid=False):
		if type == 'nxyz':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = nxyz_to_nmrdata(file)

		elif type == 'orca':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = orca_to_nmrdata(file, molid)

		elif type == 'sygcml':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = sygcml_to_nmrdata(file)

		elif type == 'jxyz':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = jxyz_to_nmrdata(file)

		elif type == 'syg_jxyz':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = syg_jxyz_to_nmrdata(file)

		elif type == 'pdb':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = pdb_to_nmrdummy(file)

		elif type == 'xyz':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.b1_coupling, self.var1b, self.b2_coupling, self.var2b, self.b3_coupling, self.var3b = xyz_to_nmrdummy(file)

	def generate_opt_in(self, prefs, path=''):

		path = path + 'optimisation/'

		self.opt_in = make_optin(prefs, str(self.molid), self.xyz, self.types, directory=path, tag='OPT')

		self.opt_log = self.opt_in.split('.')[0]
		self.opt_status = 'pre-submission'

	def generate_nmr_in(self, prefs, path=''):

		path = path + 'NMR/'

		self.nmr_in = make_nmrin(prefs, str(self.molid), self.xyz, self.types, directory=path, tag='NMR')

		self.nmr_log = self.nmr_in.split('.')[0] + '.log'
		self.nmr_status = 'pre-submission'

	def check_opt(self):
		if self.opt_status == 'None':
			print('Cannot check optimisation before in file has been created')
			return None

		self.opt_log = self.opt_log.split('.')[0]
		self.opt_status = check_opt_status(self.opt_log + '.log')

		if self.opt_status == 'successful':
			_, self.xyz, self.types, _ = generic_pybel_read(self.opt_log + '.xyz', 'xyz')
			self.energy = eread(self.opt_log + '_property.txt')

	def check_nmr(self):
		if self.nmr_status == 'None':
			print('Cannot check optimisation before in file has been created')
			return None

		self.nmr_status = check_opt_status(self.nmr_log)

		if self.nmr_status == 'Successful':
			data_from_file(self, self.nmr_log, type='g09', molid=self.molid)


	def print_xyz(self, path, num=-404):
		outname = path + str(self.molid) + '.xyz'
		nmrmol_to_xyz(self, outname, num=num)

		return outname





###

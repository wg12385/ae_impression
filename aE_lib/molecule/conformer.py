import numpy as np
import glob
from conversion.convert import nxyz_to_nmrdata, g09_to_nmrdata, optg09_to_nmrdata, sygcml_to_nmrdata, jxyz_to_nmrdata, syg_jxyz_to_nmrdata, pdb_to_nmrdummy, xyz_to_nmrdummy
from conversion.convert import nmrmol_to_xyz, nmrmol_to_nxyz
from conversion.read_g09 import eread, check_opt_status, check_nmr_status
from conversion.make_g09 import make_optcom, make_nmrcom

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
		self.opt_com = 'None'
		self.opt_log = 'None'
		self.opt_status = 'None'

		# store location and status of NMR files
		self.nmr_com = 'None'
		self.nmr_log = 'None'
		self.nmr_status = 'None'

		self.energy = 404.40404
		self.pop = 4.0404

	def get_energy_from_file(self, file):
		self.energy = eread(file)

	def data_from_file(self, file, type='g09'):
		if type == 'nxyz':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = nxyz_to_nmrdata(file)

		elif type == 'g09':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = g09_to_nmrdata(file)

		elif type == 'opt_g09':
			self.molid, self.types, self.xyz, self.dist, self.energy = opt_g09_to_nmrdata(file)

		elif type == 'sygcml':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = sygcml_to_nmrdata(file)

		elif type == 'jxyz':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = jxyz_to_nmrdata(file)

		elif type == 'syg_jxyz':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = syg_jxyz_to_nmrdata(file)

		elif type == 'pdb':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = pdb_to_nmrdummy(file)

		elif type == 'xyz':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = xyz_to_nmrdummy(file)

	def generate_opt_com(self, path='', charge=0, multiplicity=1, memory=26, processors=8,
						opt='tight', freq=True, functional='mpw1pw91', basis_set='6-31g(d)',
						solvent='none', solventmodel='iefpcm' , grid='ultrafine', direct_cmd_line_opt='none'):

		path = path + 'optimisation/'

		self.opt_com = make_optcom(str(self.molid), self.xyz, self.types, directory=path,
							charge=charge, multiplicity=multiplicity, memory=memory, processors=processors,
							opt=opt, freq=freq, functional=functional, basis_set=basis_set,
							solvent=solvent, solventmodel=solventmodel, grid=grid, direct_cmd_line_opt=direct_cmd_line_opt)

		self.opt_log = self.opt_com.split('.')[0] + '.log'
		self.opt_status = 'pre-submission'

	def generate_nmr_com(self, path='', charge=0, multiplicity=1, memory=26, processors=8,
						functional='mpw1pw91', basis_set='6-31g(d)', mixed=True,
						solvent='none', solventmodel='iefpcm', direct_cmd_line_opt='none'):

		path = path + 'NMR/'

		self.nmr_com = make_nmrcom(str(self.molid), self.xyz, self.types, directory=path,
							charge=charge, multiplicity=multiplicity, memory=memory, processors=processors,
							functional=functional, basis_set=basis_set, mixed=True,
							solvent=solvent, solventmodel=solventmodel, direct_cmd_line_opt=direct_cmd_line_opt)

		self.nmr_log = self.nmr_com.split('.')[0] + '.log'
		self.nmr_status = 'pre-submission'

	def check_opt(self):
		if self.opt_status == 'None':
			print('Cannot check optimisation before com file has been created')
			return None

		self.opt_status = check_opt_status(self.opt_log)

		if self.opt_status == 'Successful':
			data_from_file(self, self.opt_log, type='opt_g09')


	def check_nmr(self):
		if self.nmr_status == 'None':
			print('Cannot check optimisation before com file has been created')
			return None

		self.nmr_status = check_opt_status(self.nmr_log)

		if self.nmr_status == 'Successful':
			data_from_file(self, self.nmr_log, type='g09')


	def print_xyz(self, path, num=-404):
		outname = path + 'id.xyz'
		nmrmol_to_xyz(self, outname, num=num)

		return outname





###

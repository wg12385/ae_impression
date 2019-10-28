import numpy as np
import pickle
from conformational_search.conformational_search import *
from conversion.convert import *
from conversion import make_orca
from molecule.conformer import conformer as conformerclass
from .nmrfile import nmrmol
import glob

# global molecule object, contains everything for one auto-ENRICH run
class molecule(object):

	# initialise object from xyz coords and types
	def __init__(self, init_xyz=[[1,1,1],[2,2,2],[3,3,3]], init_types=[1,1,1], name='NONAME_MOL', path='', from_file=False, from_type='pickle', nosave=False):
		if from_file == False:
			# set name
			self.molid = name
			# store base xyz coordinates
			self.xyz = init_xyz
			# store types
			self.types = init_types

			# placeholder for conformer objects once defined
			self.conformers = []
			self.pops = []

			# placeholders for NMR data
			self.shift = np.zeros(len(self.types), dtype=np.float64)
			self.shift_var = np.zeros(len(self.types), dtype=np.float64)

			self.coupling1b = []
			self.var1b = []
			self.coupling2b = []
			self.var2b = []
			self.coupling3b = []
			self.var3b = []
			self.coupling4b = []
			self.var4b = []
			self.coupling5b = []
			self.var5b = []

			self.stage = 'init'
			# stages = init, pre-opt, running-opt, pre-nmr, running-nmr, post
			self.pickle = 'none'
			self.path = path

		else:
			self.pickle = 'none'
			self.read_molecule_from_file(from_file, path=path, type=from_type)

		if not nosave:
			self.save_molecule_to_file(file='', type='pickle')

	def save_molecule_to_file(self, file='', path='', type='pickle'):
		if type == 'pickle':
			if len(file) == 0:
				if len(path) == 0:
					file = self.path + self.molid + '.pkl'
				else:
					file = path + self.molid + '.pkl'
			self.pickle = file
			pickle.dump([self.pickle, self.molid, self.stage, self.xyz, self.types,
							self.conformers, self.pops, self.shift,
							self.coupling1b, self.coupling2b, self.coupling3b], open( file, "wb" ) )

		if type == 'nmredata':
			if len(file) == 0:
				file = self.path + self.molid + '.nmredata.sdf'
			nmrmol_to_nmredata(self, file)


	def read_molecule_from_file(self, file, path, type='nxyz'):
		self.stage = 'init'
		self.conformers = []
		self.pops = []
		if len(path) > 0:
			self.path = path

		filename = self.path + file

		if type == 'nxyz':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = nxyz_to_nmrdata(filename)

		if type == 'pickle':
			self.pickle, self.molid, self.stage, self.xyz, self.types, self.conformers, self.pops, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = pickle.load(open(filename, "rb" ))
			self.shift_var = np.zeros(len(self.shift), dtype=np.float64)
			self.var1b = np.zeros(len(self.coupling1b), dtype=np.float64)
			self.var2b = np.zeros(len(self.coupling2b), dtype=np.float64)
			self.var3b = np.zeros(len(self.coupling3b), dtype=np.float64)

			self.molid = self.molid.split('.')[0]

		if type == 'nmredata':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = nmredata_to_nmrmol(filename)


		elif type == 'g09':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = g09_to_nmrdata(filename)

		elif type == 'MIN_g09':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = MINIMIZE_g09_to_nmrdata(filename)

		elif type == 'tensor_g09':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = TENSOR_g09_to_nmrdata(filename)

		elif type == 'sygcml':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = sygcml_to_nmrdata(filename)

		elif type == 'pdb':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = pdb_to_nmrdummy(filename)

		elif type == 'mol2':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = mol2_to_nmrdummy(filename)

		elif type == 'xyz':
			self.molid, self.types, self.xyz, self.dist, self.shift, self.shift_var, self.coupling1b, self.var1b, self.coupling2b, self.var2b, self.coupling3b, self.var3b = xyz_to_nmrdummy(filename)


	# function to remove redundant conformers
	def remove_redundant(threshold=0.001):
		self.conformers = redundant_conformer_elim(self.conformers, threshold=threshold)

	# average NMR properties
	def boltzmann_average(self):
		self.shift = boltzmann_shift(self.conformers)
		self.coupling1b, self.coupling2b, self.coupling3b, self.coupling4b, self.coupling5b = boltzmann_coupling(self.conformers)

	# do conformational searching
	def generate_conformers(self, smiles, path='', iterations=100, RMSthresh=1, maxconfs=100, Ethresh=100000):
		xyzs, energies = torsional_search(self, smiles, path=self.path,
								iterations=iterations, RMSthresh=RMSthresh)
		print('Initial search complete,', len(xyzs), 'conformers found')

		xyzs, energies = select_conformers(xyzs, energies, maxconfs=maxconfs, Ethresh=Ethresh)

		for x, xyz in enumerate(xyzs):
			new_conf = conformerclass(xyz, self.types, x)
			new_conf.energy = energies[x]
			self.conformers.append(new_conf)

		for conformer in self.conformers:
			file = conformer.print_xyz(self.path + 'conf_search/')
			conformer.xyz_file = file


		self.stage = 'pre-opt'


	def print_xyzs(self, path=''):
		for conformer in self.conformers:
			conformer.print_xyz(path)


	# make optimsation com files
	def make_opt_in(self, prefs, path=''):

		for conformer in self.conformers:
			conformer.generate_opt_in(prefs, path=path)

		self.stage = 'running-opt'

	def make_opt_sub(self, prefs, path='', start=-1, end=-1, failed_only=False):

		innames = []
		for conformer in self.conformers:
			if conformer.opt_status != 'successful' or not failed_only:
				innames.append(conformer.opt_in)


		nodes = 1
		mem = prefs['optimisation']['memory']
		ppn = prefs['optimisation']['processors']
		walltime = prefs['optimisation']['walltime']

		in_array = make_orca.make_submission_array(self.molid, innames, path=path)
		qsub_names = make_orca.make_submission_qsub(prefs, in_array, innames, self.molid, path=path,
		 										nodes=nodes, ppn=ppn, walltime=walltime, mem=mem, start=start, end=end)

		return qsub_names

	# make nmr in files
	def make_nmr_in(self, prefs, path=''):

		for conformer in self.conformers:
			if conformer.opt_status == 'successful':
				conformer.generate_nmr_in(prefs, path=path)

		self.stage = 'running-nmr'

	def make_nmr_sub(self, prefs, path='', start=-1, end=-1, failed_only=False):

		path = path + 'NMR/'

		innames = []
		for conformer in self.conformers:
			if conformer.nmr_status != 'successful' or not failed_only:
				innames.append(conformer.nmr_in)

		nodes = 1
		mem = prefs['NMR']['memory']
		ppn = prefs['NMR']['processors']
		walltime = prefs['NMR']['walltime']

		in_array = make_orca.make_submission_array(prefs, self.molid, innames, path=path)
		qsub_names = make_orca.make_submission_qsub(prefs, in_array, innames, self.molid, path=path,
		 										nodes=nodes, ppn=ppn, walltime=walltime, mem=mem, start=start, end=end)

		return qsub_names

	# check optimistion files
	def update_opt_status(self):
		for conformer in self.conformers:
			conformer.check_opt()

	# check nmr files
	def update_nmr_status(self):
		for conformer in self.conformers:
			conformer.check_nmr()

	# compare to experimental nxyz file
	def compare_to_experiment(exp_file):
		analysis.compare_to_experiment(self, exp_file)


###

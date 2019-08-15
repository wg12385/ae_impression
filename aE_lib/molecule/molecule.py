import numpy as np
import pickle
from conformational_search.conformational_search import *
from conversion.convert import raw_to_xyz
from conversion import make_g09
from molecule.conformer import conformer as conformerclass
import glob

# global molecule object, contains everything for one auto-ENRICH run
class molecule(object):

	# initialise object from xyz coords and types
	def __init__(self, init_xyz, init_types, name='NONAME_MOL', path='', from_file=False):
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
			self.coupling1b = []
			self.coupling2b = []
			self.coupling3b = []
			self.coupling4b = []
			self.coupling5b = []

			self.stage = 'init'
			# stages = init, pre-opt, running-opt, pre-nmr, running-nmr, post
			self.path = path

			self.save_molecule(path=path)

		else:
			if from_file:
				self.molid = name
				self.load_molecule(path=path)
			else:
				self.molid = name
				self.load_molecule(path=path, overrd_filename=from_file)

	def save_molecule(self, path=''):
		filename = path + self.molid + '.pkl'
		pickle.dump([self.path, self.molid, self.stage, self.xyz, self.types,
						self.conformers, self.pops, self.shift,
						self.coupling1b, self.coupling2b, self.coupling3b], open( filename, "wb" ) )

	def load_molecule(self, path='', overrd_filename=''):
		filename = path + self.molid + '.pkl'
		if len(overrd_filename) > 1:
			filename = overrd_filename

		self.path, self.molid, self.stage, self.xyz, self.types, self.conformers, self.pops, self.shift, self.coupling1b, self.coupling2b, self.coupling3b = pickle.load(open(filename, "rb" ))


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
			file = self.path + 'conf_search/' + 'conf' + str(conformer.id) + '.xyz'
			titleline = str(conformer.id) + '  ' + str(conformer.energy)
			raw_to_xyz(file, titleline, conformer.xyz, conformer.types)
			conformer.xyz_file = file


		self.stage = 'pre-opt'

	# make optimsation com files
	def make_opt_com(self, path='', charge=0, multiplicity=1, memory=26, processors=8,
						opt='tight', freq=True, functional='mpw1pw91', basis_set='6-31g(d)',
						solvent='none', solventmodel='iefpcm' , grid='ultrafine', direct_cmd_line_opt='none'):

		for conformer in self.conformers:
			conformer.generate_opt_com(path=path, charge=charge, multiplicity=multiplicity,
										memory=memory, processors=processors, opt=opt, freq=freq,
										functional=functional, basis_set=basis_set, solvent=solvent,
										solventmodel=solventmodel, grid=grid, direct_cmd_line_opt=direct_cmd_line_opt)

		self.stage = 'running-opt'

	def make_opt_sub(self, path='', parallel=True, system='BC3', nodes=1, ppn=8,
							walltime='100:00:00', mem=26, start=-1, end=-1):

		path = path + 'optimisation/'

		comnames = []
		for conformer in self.conformers:
			comnames.append(conformer.opt_com)

		com_array = make_g09.make_submission_array(self.molid, comnames, path=path)
		qsub_names = make_g09.make_submission_qsub(com_array, comnames, self.molid, path=path, parallel=parallel, system=system,
										nodes=nodes, ppn=ppn, walltime=walltime, mem=mem, start=start, end=end)

		return qsub_names

	# make nmr com files
	def make_nmr_com(self, path='', charge=0, multiplicity=1, memory=26, processors=8,
						functional='mpw1pw91', basis_set='6-31g(d)', mixed=True,
						solvent='none', solventmodel='iefpcm', direct_cmd_line_opt='none'):

		for conformer in self.conformers:
			conformer.generate_nmr_com(path=path, charge=charge, multiplicity=multiplicity,
								memory=memory, processors=processors,
								functional=functional, basis_set=basis_set, mixed=True,
								solvent=solvent, solventmodel=solventmodel, direct_cmd_line_opt=direct_cmd_line_opt)

		self.stage = 'running-nmr'

	def make_nmr_sub(self, path='', parallel=True, system='BC3', nodes=1, ppn=8,
							walltime='100:00:00', mem=26, start=-1, end=-1):

		path = path + 'NMR/'

		comnames = []
		for conformer in self.conformers:
			comnames.append(conformer.nmr_com)

		com_array = make_g09.make_submission_array(self.molid, comnames, path=path)
		qsub_names = make_g09.make_submission_qsub(com_array, comnames, self.molid, path=path, parallel=parallel, system=system,
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

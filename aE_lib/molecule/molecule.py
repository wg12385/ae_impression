import numpy as np
import pickle
import conformational_search.conformational_search as conf_search
from molecule.conformer import conformer as conformerclass
from .nmrmol import nmrmol
import glob

# global molecule object, contains everything for one auto-ENRICH run
class molecule(nmrmol):

	def __init__(self):
		nmrmol.__init(self))

		self.conformers = []
		self.stage = 'init'

	# function to remove redundant conformers
	def remove_redundant(self, threshold=0.001):
		self.conformers = redundant_conformer_elim(self.conformers, threshold=threshold)

	# average NMR properties
	def boltzmann_average(self):
		# should there be equivalent bolzmann functions for variance ??
		self.shift = boltzmann_shift(self.conformers)
		self.coupling = boltzmann_coupling(self.conformers)
		self.variance = boltzmann_variance(self.conformers)

	# do conformational searching
	def generate_conformers(self, smiles, path='', iterations=100, RMSthresh=1, maxconfs=100, Ethresh=100000):
		xyzs, energies = conf_search.torsional_search(self, smiles, path=self.path,
								iterations=iterations, RMSthresh=RMSthresh)
		print('Initial search complete,', len(xyzs), 'conformers found')

		xyzs, energies = conf_search.select_conformers(xyzs, energies, maxconfs=maxconfs, Ethresh=Ethresh)

		for x, xyz in enumerate(xyzs):
			new_conf = conformerclass(xyz, self.types, x)
			new_conf.energy = energies[x]
			self.conformers.append(new_conf)

		for conformer in self.conformers:
			file = conformer.print_xyz(self.path + 'conf_search/')
			conformer.xyz_file = file


		self.stage = 'pre-opt'


###

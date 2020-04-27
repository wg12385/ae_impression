# Copyright 2020 Will Gerrard
#This file is part of autoenrich.

#autoenrich is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoenrich is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import pickle

from autoenrich.conformational_search import conformational_search as conf_search
from .conformer import conformer as conformerclass
from .nmrmol import nmrmol
from autoenrich.boltzmann.population import get_pop_array
from autoenrich.boltzmann.averaging import *
import glob
from tqdm import tqdm

# global molecule object, contains everything for one auto-ENRICH run
class molecule(nmrmol):

	def __init__(self, molid, path=''):
		nmrmol.__init__(self, molid, path)

		self.conformers = []
		self.stage = 'init'

	# average NMR properties
	def boltzmann_average(self):

		pops = get_pop_array(self.conformers)
		for c in range(len(self.conformers)):
			self.conformers[c].pop = pops[c]

		self.shift, self.shift_var = boltzmann_shift(self.conformers)
		self.coupling, self.coupling_var = boltzmann_coupling(self.conformers)

	# do conformational searching
	def generate_conformers(self, smiles, path='', iterations=100, RMSthresh=1, maxconfs=100, Ethresh=100000):
		xyzs, energies = conf_search.torsional_search(smiles, iterations=iterations, RMSthresh=RMSthresh)
		print('Initial search complete,', len(xyzs), 'conformers found')

		xyzs, energies = conf_search.select_conformers(xyzs, energies, maxconfs=maxconfs, Ethresh=Ethresh)

		for x, xyz in enumerate(xyzs):
			new_conf = conformerclass(str(x), path=path)
			new_conf.xyz = xyz
			new_conf.types = self.types
			new_conf.conn = self.conn
			new_conf.coupling_len = self.coupling_len
			new_conf.energy = energies[x]
			self.conformers.append(new_conf)

		self.stage = 'pre-opt'

	def remove_redundant(self, RMS_thresh=3.0, MMe_thresh=100000, DFTe_thresh=100000, maxconfs=200):
		size = len(self.conformers[0].types)
		confs = len(self.conformers)
		redundant = []

		print('Calculating distance matrices')
		for conf in tqdm(self.conformers):
			conf.get_distance_matrix()
			conf.redundant = False

		print('Removing conformers based on energy threshold')
		for conf in tqdm(self.conformers):
			if conf.redundant:
				continue

			if conf.opt_status == 'successful':
				if conf.energy > DFTe_thresh:
					redundant.append(conf.molid)
					continue
			else:
				if conf.energy > MMe_thresh:
					redundant.append(conf.molid)
					continue

		print('Removing conformers based on RMS threshold')
		dist = np.zeros((confs, confs), dtype=np.float64)
		c1 = -1
		for conf1 in tqdm(self.conformers):
			c1 += 1
			c2 = -1
			for conf2 in self.conformers:
				c2 += 1
				if c1 >= c2:
					continue

				if conf2.molid not in redundant:
					dist[c1][c2] = np.sqrt(np.mean(np.square(conf1.dist-conf2.dist)))
					if dist[c1][c2] > RMS_thresh:
						redundant.append(conf2.molid)


		if (len(self.conformers) - len(set(redundant))) > maxconfs:
			print('Still too many conformers, selecting least similar')
			to_remove = len(self.conformers) - maxconfs
			# Keep looping until enough conformers are marked for removal
			with tqdm(total=to_remove-len(set(redundant))) as pbar:
				while len(redundant) < to_remove:
					id = 0
					lowest_dist = 9999999999999999
					# Loop over conformers
					for c, conf in enumerate(self.conformers):
						# Get sum of distances to all other conformers for conformer i
						sumdist = np.sum(dist[c])
						# Store lowest distance and conformer id
						if sumdist <= lowest_dist and conf.molid not in redundant:
							id = conf.molid
							lowest_dist = sumdist
					# add least geometrically different conformer to removal list
					redundant.append(id)
					pbar.update(1)

		conformers = len(self.conformers)
		for conf in self.conformers:
			if conf.molid in redundant:
				conf.redundant = True
				conformers -= 1

		print('Conformers removed, number of conformers to be used: ', conformers)














###

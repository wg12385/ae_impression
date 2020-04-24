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

	# function to remove redundant conformers
	def remove_redundant(self, threshold=0.001):
		self.conformers = redundant_conformer_elim(self.conformers, threshold=threshold)

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

	def remove_redundant(self, RMS_thresh=3.0, MMe_thresh=100000, DFTe_thresh=100000, maxconfs=1000):
		size = len(self.conformers[0].types)

		print('Making distance matrices')
		for conf in tqdm(self.conformers):
			dist1 = np.zeros((size, size), dtype=np.float64)
			for i in np.where(conf.types != 1)[0]:
				for j in np.where(conf.types != 1)[0]:
					dist1[i][j] = np.sqrt(np.square(conf.xyz[i][0]-conf.xyz[j][0])
										+ np.square(conf.xyz[i][1]-conf.xyz[j][1])
										+ np.square(conf.xyz[i][2]-conf.xyz[j][2]))
			conf.dist = dist1

		for RMS_thresh in range(1, 40):
			RMS_thresh = RMS_thresh * 0.1

			redundant = []
			c1 = 0
			for conf1 in tqdm(self.conformers):
				c1 += 1
				conf1.redundant = False

				if conf1.opt_status == 'successful':
					if conf1.energy > DFTe_thresh:
						redundant.append(conf1.molid)
						continue
				else:
					if conf1.energy > MMe_thresh:
						redundant.append(conf1.molid)
						continue

				for c2, conf2 in enumerate(self.conformers):
					if c1 >= c2:
						continue
					if conf2.molid in redundant:
						continue

					RMS = np.sqrt(np.mean(np.square(conf1.dist-conf2.dist)))
					#RMS = np.sqrt(np.mean(np.square(np.asarray(conf1.xyz)-np.asarray(conf2.xyz))))

					if RMS > RMS_thresh:
						redundant.append(conf2.molid)

			'''
			for conf in self.conformers:
				if conf.molid in redundant:
					conf.redundant = True
			'''
			print('RMS: ', '{0:<3.1f}'.format(RMS_thresh) ,' number of conformers:', len(self.conformers) - len(set(redundant)))
			'''
			red = 0
			for conf in self.conformers:
				if conf.redundant:
					red += 1

			print('Redundant: ', red)
			'''














###

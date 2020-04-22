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


from .nmrmol import nmrmol
from autoenrich.util.file_gettype import get_type
from autoenrich.util.flag_handler.hdl_targetflag import flag_to_target
from autoenrich.util.filename_utils import get_unique_part
from autoenrich.ml.features import GNR_features

from tqdm import tqdm
import gzip
import pickle

import numpy as np
np.set_printoptions(threshold=99999999999)

class dataset(object):

	def __init__(self):

		self.mols = []
		self.x = []
		self.y = []
		self.r = []
		self.mol_order = []
		self.big_data = False
		self.files = []

		self.params = {}


	def get_mols(self, files, type='', label_part=-1, fallback=False):
		self.mols = []
		self.files = files

		if len(files) > 2000:
			self.big_data = True

		if label_part == -1:
			label_part = get_unique_part(files)
			if label_part == -1:
				fallback = True

		for f, file in enumerate(files):
			if fallback:
				id = str(f)
			else:
				id = file.split('/')[-1].split('.')[0].split('_')[label_part]

			try:
				int(id)
				id = 'nmrmol' + str(id)
			except:
				pass

			if self.big_data:
				self.mols.append([file, id, type])
			else:
				mol = nmrmol(molid=id)

				if type == '':
					ftype = get_type(file)
				else:
					ftype = type
				mol.read_nmr(file, ftype)
				self.mols.append(mol)


	def get_features_frommols(self, args, params={}, molcheck_run=False, training=True):

		featureflag = args['featureflag']
		targetflag = args['targetflag']
		try:
			max = args['max_size']
		except:
			max = 200

		for mol in self.mols:
			if len(mol.types) > max:
				max = len(mol.types)
				print('WARNING, SETTING MAXIMUM MOLECULE SIZE TO, ', max)

		if 'cutoff' in params:
			if params['cutoff'] < 0.1:
				params['cutoff'] = 0.1
		else:
			params['cutoff'] = 5.0

		x = []
		y = []
		r = []

		self.params = params

		target = flag_to_target(targetflag)
		self.remove_mols(target)
		if molcheck_run:
			return 

		if featureflag in ['aSLATM', 'CMAT', 'FCHL', 'ACSF']:
			import qml
		elif featureflag in ['BCAI']:
			from autoenrich.ml.features import TFM_features

		_, y, r = GNR_features.get_dummy_features(self.mols, targetflag)

		if featureflag == 'aSLATM':
			mbtypes = [[1],[1,1], [1,1,1], [1,1,6], [1,1,7], [1,1,8], [1,1,9], [1,6], [1,6,1], [1,6,6], [1,6,7], [1,6,8], [1,6,9], [1,7], [1,7,1], [1,7,6], [1,7,7], [1,7,8], [1,7,9], [1,8], [1,8,1], [1,8,6], [1,8,7], [1,8,8], [1,8,9], [1,9], [1,9,1], [1,9,6], [1,9,7], [1,9,8], [1,9,9], [6], [6,1], [6,1,1], [6,1,6], [6,1,7], [6,1,8], [6,1,9], [6,6], [6,6,1], [6,6,6], [6,6,7], [6,6,8], [6,6,9], [6,7], [6,7,1], [6,7,6], [6,7,7], [6,7,8], [6,7,9], [6,8], [6,8,1], [6,8,6], [6,8,7], [6,8,8], [6,8,9], [6,9], [6,9,1], [6,9,6], [6,9,7], [6,9,8], [6,9,9], [7],[7,1], [7,1,1], [7,1,6], [7,1,7], [7,1,8], [7,1,9], [7,6], [7,6,1], [7,6,6], [7,6,7], [7,6,8], [7,6,9], [7,7], [7,7,1], [7,7,6], [7,7,7], [7,7,8], [7,7,9], [7,8], [7,8,1], [7,8,6], [7,8,7], [7,8,8], [7,8,9], [7,9], [7,9,1], [7,9,6], [7,9,7], [7,9,8], [7,9,9], [8], [8,1], [8,1,1], [8,1,6], [8,1,7], [8,1,8], [8,1,9], [8,6], [8,6,1], [8,6,6], [8,6,7], [8,6,8], [8,6,9], [8,7], [8,7,1], [8,7,6], [8,7,7], [8,7,8], [8,7,9], [8,8], [8,8,1], [8,8,6], [8,8,7], [8,8,8], [8,8,9], [8,9], [8,9,1], [8,9,6], [8,9,7], [8,9,8], [8,9,9], [9], [9,1], [9,1,1], [9,1,6], [9,1,7], [9,1,8], [9,1,9], [9,6], [9,6,1], [9,6,6], [9,6,7], [9,6,8], [9,6,9], [9,7], [9,7,1], [9,7,6], [9,7,7], [9,7,8], [9,7,9], [9,8], [9,8,1], [9,8,6], [9,8,7], [9,8,8], [9,8,9], [9,9], [9,9,1], [9,9,6], [9,9,7], [9,9,8], [9,9,9]]
			'''
			nuclear_charges = []
			for tmp_mol in mols:
				nuclear_charges.append(tmp_mol.types)
			mbtypes = qml.representations.get_slatm_mbtypes(nuclear_charges)
			'''
			reps = qml.representations.generate_slatm(mol.xyz, mol.types, mbtypes, rcut=cutoff)
			x = np.asarray(reps)

		elif featureflag == 'CMAT':
			reps = qml.representations.generate_atomic_coulomb_matrix(mol.types, mol.xyz, size = max, central_cutoff = cutoff)
			x = np.asarray(reps)

		elif featureflag == 'FCHL':
			reps = qml.fchl.generate_representation(mol.xyz, mol.types, max, cut_distance=cutoff)
			x = np.asarray(reps)

		elif featureflag == 'ACSF':
			reps = qml.representations.generate_acsf(mol.types, mol.xyz, elements=[1, 6, 7, 8, 9, 14, 15, 16, 17, 35],
													nRs2=int(nRs2), nRs3=int(nRs3),
													nTs=int(nTs), eta2=eta2, eta3=eta3, zeta=zeta, rcut=cutoff, acut=acut,
													bin_min=0.0, gradients=False)
			x = np.asarray(reps)

		elif featureflag == 'BCAI':

				_x, _y, _r, mol_order = TFM_features.get_BCAI_features(self.mols, targetflag, training=training)

				x.extend(_x)
				y.extend(_y)
				r.extend(_r)
				batch_mols = []

		else:
			print('Feature flag not recognised, no feature flag: ', featureflag)

		if featureflag == 'BCAI':
			self.x = x
			self.y = y
			self.r = r
			self.mol_order = mol_order
		else:
			self.x = np.asarray(x)
			self.y = np.asarray(y)
			self.r = r

		if featureflag not in ['dummy', 'BCAI']:
			print('Reps generated, shape: ', self.x.shape)



	def assign_from_ml(self, pred_y, var, zero=True):
		assert len(self.r) > 0, print('Something went wrong, nothing to assign')

		try:
			with gzip.open(self.r[0], "rb") as f:
				self.r = pickle.load(f)
		except Exception as e:
			print(e)
			pass

		for molrf in tqdm(self.mols, desc='Assigning predictions'):

			if type(molrf) == list:
					mol = nmrmol(molid=molrf[1])

					if molrf[2] == '':
						ftype = get_type(molrf[2])
					else:
						ftype = molrf[2]

					mol.read_nmr(molrf[0], ftype)
			else:
				mol = molrf

			if zero:
				mol.coupling = np.zeros((len(mol.types), len(mol.types)), dtype=np.float64)
				mol.shift = np.zeros((len(mol.types)), dtype=np.float64)

			for m, refs in enumerate(self.r):
				for r, ref in enumerate(refs):
					if mol.molid != ref[0]:
						continue

					if len(ref) == 2:
						for t in range(len(mol.types)):
							if ref[1] == t:
								mol.shift[t] = pred_y[r]
								mol.shift_var[t] = var[r]

					elif len(ref) == 3:
						for t1 in range(len(mol.types)):
							for t2 in range(len(mol.types)):
								if ref[1] == t1 and ref[2] == t2:
									mol.coupling[t1][t2] = pred_y[r]
									mol.coupling_var[t1][t2] = var[r]

	def remove_mols(self, target):
		## Discard useless molecules:
		to_remove = []
		print('Checking structures')
		for molrf in tqdm(self.mols):
			if self.big_data:
				mol = nmrmol(molid=molrf[1])

				if molrf[2] == '':
					ftype = get_type(molrf[0])
				else:
					ftype = molrf[2]

				mol.read_nmr(molrf[0], ftype)
			else:
				mol = molrf

			if len(target) == 1:
				if not target in mol.types:
					to_remove.append(mol.molid)
			elif len(target) == 3:
				found = False
				for i, it in enumerate(mol.types):
					for j, jt in enumerate(mol.types):
						if i >= j:
							continue

						if mol.coupling_len[i][j] == target[0]:
							if it == target[1] and jt == target[2]:
								found = True
							elif it == target[2] and jt == target[1]:
								found = True

				if not found:
					to_remove.append(mol.molid)

		if len(to_remove) > 1:
			keep = []
			for i in range(len(self.mols)):
				if self.big_data:
					if self.mols[i][1] in to_remove:
						continue
				else:
					if self.mols[i].molid in to_remove:
						continue
				keep.append(self.mols[i])

			self.mols = keep

			print('REMOVED ', len(to_remove), '/', len(to_remove)+len(keep), ' molecules due to lack of features')
			#print(to_remove[:10])
			assert len(keep) > 1

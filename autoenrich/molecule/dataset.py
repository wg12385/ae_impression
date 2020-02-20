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

from tqdm import tqdm

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


	def get_features_frommols(self, args, params={}, molcheck_run=False):

		featureflag = args['featureflag']
		targetflag = args['targetflag']
		try:
			max = args['max_size']
		except:
			max = 200

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
			print('molcheck complete')
			return

		if featureflag in ['aSLATM', 'CMAT', 'FCHL', 'ACSF']:
			from autoenrich.ml.features import QML_features
		elif featureflag in ['BCAI']:
			from autoenrich.ml.features import TFM_features
		elif featureflag in ['dummy']:
			from autoenrich.ml.features import GNR_features

		if featureflag == 'dummy':
			for mol in self.mols:
				_x, _y, _r = GNR_features.get_dummy_features([mol], targetflag)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)


		elif featureflag == 'aSLATM':
			mbtypes = QML_features.get_aSLATM_mbtypes(self.mols)
			for mol in self.mols:
				_x, _y, _r = QML_features.get_aSLATM_features([mol], targetflag, params['cutoff'], max=max, mbtypes=mbtypes)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)


		elif featureflag == 'CMAT':

			# Set (not found) parameters to defaults
			if not 'central_decay' in params:
				params['central_decay'] = -1
			if not 'interaction_cutoff' in params:
				params['interaction_cutoff'] = 1000000.0
			if not 'interaction_decay' in params:
				params['interaction_decay'] = -1

			for mol in self.mols:
				if len(mol.types) >max:
					args['max_size'] = len(mol.types)
					print('WARNING, SETTING MAXIMUM MOLECULE SIZE TO, ', max)
			for mol in self.mols:
				_x, _y, _r = QML_features.get_CMAT_features([mol], targetflag, params['cutoff'], args['max_size'], central_decay=params['central_decay'],
														interaction_cutoff=params['interaction_cutoff'], interaction_decay=params['interaction_decay'])
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)


		elif featureflag == 'FCHL':
			for mol in self.mols:
				if len(mol.types) >max:
					max = len(mol.types)
					print('WARNING, SETTING MAXIMUM MOLECULE SIZE TO, ', max)

				'''
				for mol in self.mols:
				_x, _y, _r = QML_features.get_FCHL_features([mol], targetflag, params['cutoff'], max)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)
				'''
			x, y, r = QML_features.get_FCHL_features(self.mols, targetflag, params['cutoff'], max)


		elif featureflag == 'ACSF':
			elements = set()
			for tmp_mol in self.mols:
				elements = elements.union(tmp_mol.types)
			elements = sorted(list(elements))

			# Set (not found) parameters to defaults


			for mol in self.mols:
				_x, _y, _r = QML_features.get_ACSF_features([mol], targetflag, params['cutoff'], elements=elements, nRs2=params['nRs2'],
												nRs3=params['nRs3'], nTs=params['nTs'], eta2=params['eta2'], eta3=params['eta3'], zeta=params['zeta'],
												acut=params['acut'], bin_min=params['bin_min'])
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)

		elif featureflag == 'BCAI':

				_x, _y, _r, mol_order = TFM_features.get_BCAI_features(self.mols, targetflag)

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

		for mol in self.mols:
			if zero:
				mol.coupling = np.zeros((len(mol.types), len(mol.types)), dtype=np.float64)
				mol.shift = np.zeros((len(mol.types)), dtype=np.float64)

			for r, ref in enumerate(self.r):

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









##

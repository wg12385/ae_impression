
from .nmrmol import nmrmol
from ml.features import QML_features, TFM_features
import numpy as np
np.set_printoptions(threshold=99999999999)

class dataset(object):

	def __init__(self):

		self.mols = []
		self.x = []
		self.y = []
		self.r = []

		self.params = {}


	def get_mols(self, files, type):
		self.mols = []
		for file in files:
			id = file.split('/')[-1].split('_')[0]
			mol = nmrmol(molid=id)
			mol.read_nmr(file, type)

			self.mols.append(mol)
			#print(id)
			#print(mol.coupling_len)


	def get_features_frommols(self, args, params={}):

		featureflag = args['featureflag']
		targetflag = args['targetflag']
		try:
			max = args['max_size']
		except:
			max = 200

		x = []
		y = []
		r = []

		self.params = params
		if featureflag == 'aSLATM':
			mbtypes = QML_features.get_aSLATM_mbtypes(self.mols)
			for mol in self.mols:
				_x, _y, _r = QML_features.get_aSLATM_features([mol], targetflag, params['cutoff'], max=max, mbtypes=mbtypes)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)


		elif featureflag == 'CMAT':
			for mol in self.mols:
				if len(mol.types) >max:
					args['max_size'] = len(mol.types)
					print('WARNING, SETTING MAXIMUM MOLECULE SIZE TO, ', max)
			for mol in self.mols:
				_x, _y, _r = QML_features.get_CMAT_features([mol], targetflag, params['cutoff'],args['max_size'], central_decay=params['central_decay'],
														interaction_cutoff=params['interaction_cutoff'], interaction_decay=params['interaction_decay'])
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)


		elif featureflag == 'FCHL':
			for mol in self.mols:
				if len(mol.types) >max:
					max = len(mol.types)
					print('WARNING, SETTING MAXIMUM MOLECULE SIZE TO, ', max)
			for mol in self.mols:
				_x, _y, _r = QML_features.get_FCHL_features([mol], targetflag, params['cutoff'], max)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)


		elif featureflag == 'ACSF':
			elements = set()
			for tmp_mol in self.mols:
				elements = elements.union(tmp_mol.types)
			elements = sorted(list(elements))
			for mol in self.mols:
				_x, _y, _r = QML_features.get_ACSF_features([mol], targetflag, params['cutoff'], elements=elements, nRs2=params['nRs2'],
												nRs3=params['nRs3'], nTs=params['nTs'], eta2=params['eta2'], eta3=params['eta3'], zeta=params['zeta'],
												acut=params['acut'], bin_min=params['bin_min'])
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)

		elif featureflag == 'BCAI':
			x, y, r = TFM_features.get_BCAI_features(self.mols, targetflag)

		else:
			print('Feature flag not recognised, no feature flag: ', featureflag)

		self.x = x
		self.y = np.asarray(y)
		self.r = r

	def get_features_fromfiles(self, files, featureflag='CMAT', targetflag='CCS', cutoff=5.0, max=400, mbtypes=[], elements=[]):
		self.params = {}

		x = []
		y = []
		r = []

		if featureflag == 'aSLATM':
			mbtypes = QML_features.get_aSLATM_mbtypes([])
			for file in files:
				id = file.split('.')[0].split('_')[-1]
				mol = nmrmol(id)
				mol.read_nmr(file, 'nmredata')
				_x, _y, _r = QML_features.get_aSLATM_features(mol, targetflag, cutoff, mbtypes)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)


		elif featureflag == 'CMAT':
			for file in files:
				id = file.split('.')[0].split('_')[-1]
				mol = nmrmol(id)
				mol.read_nmr(file, 'nmredata')
				_x, _y, _r = QML_features.get_CMAT_features(mol, targetflag, cutoff, max)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)


		elif featureflag == 'FCHL':
			for file in files:
				id = file.split('.')[0].split('_')[-1]
				mol = nmrmol(id)
				mol.read_nmr(file, 'nmredata')
				_x, _y, _r = QML_features.get_FCHL_features(mol, targetflag, cutoff, max)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)


		elif featureflag == 'ACSF':
			for file in files:
				id = file.split('.')[0].split('_')[-1]
				mol = nmrmol(id)
				mol.read_nmr(file, 'nmredata')
				_x, _y, _r = QML_features.get_ACSF_features(mol, targetflag, cutoff, elements)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)

		elif featureflag == 'BCAI':
			for file in files:
				id = file.split('.')[0].split('_')[-1]
				mol = nmrmol(id)
				mol.read_nmr(file, 'nmredata')
				_x, _y, _r = TFM_features.get_BCAI_features(mol, targetflag, cutoff)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)

		self.x = np.asarray(x)
		self.y = np.asarray(y)
		self.r = r


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











##

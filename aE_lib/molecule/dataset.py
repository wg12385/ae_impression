
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
			id = file.split('.')[0].split('_')[-1]
			mol = nmrmol(molid=id)
			mol.read_nmr(file, type)
			self.mols.append(mol)
			#print(id)
			#print(mol.coupling_len)


	def get_features_frommols(self, featureflag='CMAT', targetflag='CCS', max=400, params={}):

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
				if len(mol.types) > max:
					max = len(mol.types)
			for mol in self.mols:
				_x, _y, _r = QML_features.get_CMAT_features([mol], targetflag, params['cutoff'], max)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)


		elif featureflag == 'FCHL':
			for mol in self.mols:
				if len(mol.types) > max:
					max = len(mol.types)
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
				_x, _y, _r = QML_features.get_ACSF_features([mol], targetflag, params['cutoff'], elements=elements)
				x.extend(_x)
				y.extend(_y)
				r.extend(_r)

		elif featureflag == 'BCAI':
			x, y, r = TFM_features.get_BCAI_features(mol, targetflag, params['cutoff'])

		else:
			print('Feature flag not recognised, no feature flag: ', featureflag)

		self.x = np.asarray(x)
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














##

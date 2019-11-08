
from .nmrmol import nmrmol
from ml.features import QML_features, TFM_features

class dataset(object):

	def __init__(self):

		self.mols = []
		self.x = []
		self.y = []
		self.r = []


	def get_mols(self, files, type):
		self.mols = []
		for file in files:
			id = file.split('.')[0].split('_')[-1]
			mol = nmrmol(id=id)
			mol.read_nmr(file)
			self.mols.append(mol)


	def get_features_frommols(self, featureflag='CMAT', targetflag='CCS', cutoff=5.0):

		if featureflag == 'aSLATM':
			mbtypes = QML_features.get_aSLATM_mbtypes(mols)
			for mol in self.mols:
				x, y, r = QML_features.get_aSLATM_features(mol, targetflag, cutoff, mbtypes)
				x.extend(x)
				y.extend(y)
				r.extend(r)


		elif featureflag == 'CMAT':
			for mol in self.mols:
				if len(mol.types) > max:
					max = len(mol.types)
			for mol in self.mols:
				x, y, r = QML_features.get_CMAT_features(mol, targetflag, cutoff, max)
				x.extend(x)
				y.extend(y)
				r.extend(r)


		elif featureflag == 'FCHL':
			for mol in self.mols:
				if len(mol.types) > max:
					max = len(mol.types)
			for mol in self.mols:
				x, y, r = QML_features.get_FCHL_features(mol, targetflag, cutoff, max)
				x.extend(x)
				y.extend(y)
				r.extend(r)


		elif featureflag == 'ACSF':
			elements = set()
			for tmp_mol in self.mols:
				elements = elements.union(tmp_mol.types)
			elements = sorted(list(elements))
			for mol in self.mols:
				x, y, r = QML_features.get_ACSF_features(mol, targetflag, cutoff)
				x.extend(x)
				y.extend(y)
				r.extend(r)

		elif featureflag == 'BCAI':
			x, y, r = TFM_features.get_BCAI_features(mol, targetflag, cutoff)

		self.x = np.asarray(x)
		self.y = np.asarray(y)
		self.r = r

	def get_features_fromfiles(self, files, featureflag='CMAT', targetflag='CCS', cutoff=5.0, max=400, mbtypes=[], elements=[]):

		x = []
		y = []
		r = []

		if featureflag == 'aSLATM':
			mbtypes = QML_features.get_aSLATM_mbtypes([])
			for file in files:
				id = file.split('.')[0].split('_')[-1]
				mol = nmrmol(id)
				mol.read_nmr(file, 'nmredata')
				x, y, r = QML_features.get_aSLATM_features(mol, targetflag, cutoff, mbtypes)
				x.extend(x)
				y.extend(y)
				r.extend(r)


		elif featureflag == 'CMAT':
			for file in files:
				id = file.split('.')[0].split('_')[-1]
				mol = nmrmol(id)
				mol.read_nmr(file, 'nmredata')
				x, y, r = QML_features.get_CMAT_features(mol, targetflag, cutoff, max)
				x.extend(x)
				y.extend(y)
				r.extend(r)


		elif featureflag == 'FCHL':
			for file in files:
				id = file.split('.')[0].split('_')[-1]
				mol = nmrmol(id)
				mol.read_nmr(file, 'nmredata')
				x, y, r = QML_features.get_FCHL_features(mol, targetflag, cutoff, max)
				x.extend(x)
				y.extend(y)
				r.extend(r)


		elif featureflag == 'ACSF':
			for file in files:
				id = file.split('.')[0].split('_')[-1]
				mol = nmrmol(id)
				mol.read_nmr(file, 'nmredata')
				x, y, r = QML_features.get_ACSF_features(mol, targetflag, cutoff, elements)
				x.extend(x)
				y.extend(y)
				r.extend(r)

		elif featureflag == 'BCAI':
			for file in files:
				id = file.split('.')[0].split('_')[-1]
				mol = nmrmol(id)
				mol.read_nmr(file, 'nmredata')
				x, y, r = TFM_features.get_BCAI_features(mol, targetflag, cutoff)
				x.extend(x)
				y.extend(y)
				r.extend(r)

		self.x = np.asarray(x)
		self.y = np.asarray(y)
		self.r = r














##

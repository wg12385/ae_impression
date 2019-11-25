import pickle
from file_creation.HPC_submission import make_HPC_header
from molecule.dataset import dataset
from file_creation.structure_formats.nmredata import nmrmol_to_nmredata
import sys
import numpy as np
import glob

def setup_predict(args):

	# NEED TO MAKE GPU VERSION
	python_env = args['python_env']
	system = args['system']
	nodes = 1
	ppn = args['processors']
	walltime = args['walltime']
	mem = args['memory']

	jobname = 'IMP_' + 'predict'
	submission_file = jobname + '.submit'
	header = make_HPC_header(jobname, system=args['system'], nodes=1, ppn=args['processors'], walltime=args['walltime'], mem=args['memory'])
	strings = []
	strings.append('')
	strings.append('shopt -s expand_aliases')
	strings.append('source ~/.bashrc')
	strings.append("source activate {env:<10s}".format(env=python_env))
	strings.append("IMPRESSION predict --prefs {0:<10s}".format(args['prefs']))


	with open(submission_file, 'w') as f:
		for string in header:
			print(string, file=f)
		for string in strings:
			print(string, file=f)


def predict(args):

	for files_set in args['test_sets']:
		parts = files_set.split('/')
		path = ''
		for part in parts[:-1]:
			path = path + part + '/'

		files = glob.glob(files_set)
		if len(files) == 0:
			print ('No file(s) found matching ', args['training_set'])
			sys.exit(0)
		dset = dataset()
		dset.get_mols(files, type='nmredata')
		if len(dset.mols) == 0:
			print('No molecules loaded. . .')
			sys.exit(0)

		for m, model_file in enumerate(args['models']):

			print('Predicting from model: ', model_file)

			model = pickle.load(open(model_file, 'rb'))

			dset.get_features_frommols(featureflag=model.args['featureflag'],
						targetflag=model.args['targetflag'], params=model.params)
			assert len(dset.x) > 0, print('No features made. . . ')

			if args['store_datasets']:
				pickle.dump(dset, open('OPT_training_set.pkl', 'wb'))

			y_pred = model.predict(dset.x)
			assert len(y_pred) == len(dset.y)

			v_preds = []
			for i in range(args['var']):
				var_model_file = model_file.split('.pkl')[0] + '_' + str(i+1) + '.pkl'

				try:
					var_model = pickle.load(open(var_model_file, 'rb'))
				except Exception as e:
					print(e)
					continue

				assert model.args['featureflag'] == var_model.args['featureflag']
				assert model.args['targetflag'] == var_model.args['targetflag']
				assert model.params == var_model.params

				print('\tPredicting from ', var_model_file)
				tmp_preds = var_model.predict(dset.x)
				v_preds.append(tmp_preds)

			if args['var'] > 0:
				var = np.var(np.asarray(v_preds), axis=0)
			else:
				var = np.zeros(len(y_pred), dtype=np.float64)


			## Currently not being zero'd properly, and values look really wrong.

			if m == 0:
				dset.assign_from_ml(y_pred, var, zero=True)
			else:
				dset.assign_from_ml(y_pred, var, zero=False)

		for mol in dset.mols:
			outname = path + 'IMP_' + mol.molid + '.nmredata.sdf'
			nmrmol_to_nmredata(mol, outname)

	print('Done')



















#####

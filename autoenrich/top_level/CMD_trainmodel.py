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


from autoenrich.file_creation.HPC_submission import make_HPC_header

import pickle
import glob
import sys

# creates a submission file that will run a hyper parameter search
def setup_trainmodel(args):

	# NEED TO MAKE GPU VERSION
	python_env = args['python_env']
	system = args['system']
	nodes = 1
	ppn = args['processors']
	walltime = args['walltime']
	mem = args['memory']

	jobname = 'IMP_' + args['modelflag'] + '_' + args['featureflag'] + '_' + args['targetflag'] + '_' + args['searchflag'] + '_HPS'
	submission_file = jobname + '.submit'
	header = make_HPC_header(jobname, system=args['system'], nodes=1, ppn=args['processors'], walltime=args['walltime'], mem=args['memory'])
	strings = []
	strings.append('')
	strings.append('shopt -s expand_aliases')
	strings.append('source ~/.bashrc')
	strings.append("source activate {env:<10s}".format(env=python_env))
	strings.append("IMPRESSION train --prefs {0:<10s}".format(args['prefs']))


	with open(submission_file, 'w') as f:
		for string in header:
			print(string, file=f)
		for string in strings:
			print(string, file=f)


def trainmodel(args):

	from autoenrich.ml.HPS import HPS
	from autoenrich.molecule.dataset import dataset
	from autoenrich.molecule.nmrmol import nmrmol

	params = []
	for param in args['param_ranges'].keys():
		params.append(param)
	args['param_list'] = params

	args['BCAI_load'] = 'True'

	if args['featureflag'] == 'BCAI':

		if args['training_set'].split('.')[-1] == 'pkl':
			dset = pickle.load(open(args['training_set'], 'rb'))
		else:
			files = glob.glob(args['training_set'])
			dset = dataset()
			dset.get_mols(files, type='nmredata')

		print('Number of molecules in training set: ', len(dset.mols))

		if args['BCAI_load'] != 'True' or len(dset.x) == 0:
			dset.get_features_frommols(args, params={}, training=True)

		with open('training_data/BCAI_dataset.pkl', "wb") as f:
			pickle.dump(dset, f)

		dset, score = HPS(dset, args)

		print('Final Score: ', score)

	else:
		if args['feature_optimisation'] == 'True':
			if args['training_set'].split('.')[-1] == 'pkl':
				dset = pickle.load(open(args['training_set'], 'rb'))
			elif args['training_set'].split('.')[-1] == 'csv':
				dset = load_dataset_from_csv(args['training_set'])
			else:
				args['load_dataset'] = 'false'
				if args['load_dataset'] == 'true':
					dset = pickle.load(open('training_data/empty_dataset.pkl', 'rb'))
				else:
					files = glob.glob(args['training_set'])
					dset = dataset()
					dset.get_mols(files, type='nmredata')
					assert len(dset.mols) > 0


			dset, score = HPS(dset, args)

			if args['store_datasets'] == 'True':
				pickle.dump(dset, open('OPT_training_set.pkl', 'wb'))


		else:
			dset = pickle.load(open(args['feature_file'], "rb"))
			assert len(dset.x) > 0
			assert len(dset.y) > 0

		HPS(dset, args)





###

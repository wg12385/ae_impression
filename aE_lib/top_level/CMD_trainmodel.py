
from file_creation.HPC_submission import make_HPC_header
from ml.HPS import HPS
from molecule.dataset import dataset
from molecule.nmrmol import nmrmol
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

	if args['feature_optimisation'] == 'True':
		if args['training_set'].split('.')[-1] == 'pkl':
			dset = pickle.load(open(args['training_set'], 'rb'))
		elif args['training_set'].split('.')[-1] == 'csv':
			dset = load_dataset_from_csv(args['training_set'])
		else:
			files = glob.glob(args['training_set'])
			dset = dataset()
			dset.get_mols(files, type='nmredata')
			assert len(dset.mols) > 0

			if args['store_datasets'] == 'True':
				pickle.dump(dset, open('training_set.pkl', 'wb'))

		dset, score = HPS(dset, args)

		if args['store_datasets']:
			pickle.dump(dset, open('OPT_training_set.pkl', 'wb'))


	else:
		dset = pickle.load(open(args['feature_file'], "rb"))
		assert len(dset.x) > 0
		assert len(dset.y) > 0

		HPS(dset, args)





###

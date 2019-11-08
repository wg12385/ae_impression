# master script to do everything IMPRESSION-y

import argparse
import sys
sys.path.append('aE_lib/')

from util.flag_handler import targetflag, flag_combos

from top_level import CMD_trainmodel, CMD_predict, CMD_compare

from util.header import print_header_IMP



if __name__ == "__main__":
	# Argparser
	parser = argparse.ArgumentParser(description='IMPRESSION')
	# Define command to perform
	parser.add_argument('Command', help='IMPRESSION command',
							choices=['setup_train', 'train', 'setup_predict', 'predict', 'compare', 'test'])

	parser.add_argument('--training_set', help='Training dataset file(s), either single csv/pkl file or one/multiple nmredata files',
	 					action="store", dest='training_set', default='None')
	parser.add_argument('--test_set1', help='Testing dataset 1',
	 					action="store", dest='test_set1', default='None')
	parser.add_argument('--test_set2', help='Testing dataset 2',
	 					action="store", dest='test_set2', default='None')
	parser.add_argument('--test_set3', help='Testing dataset 3',
	 					action="store", dest='test_set3', default='None')

	parser.add_argument('--store_datasets', help='Option to store datasets as pickle files for later use',
	 					action="store", dest='store_datasets', default='False')

	# Optional arguments for train/test
	parser.add_argument('--modelflag', help='type of model to use',
							action="store", dest='modelflag', default='KRR',
							choices=['KRR', 'FCHL', 'NN', 'TFM'])
	parser.add_argument('--featureflag', help='type of features to use',
							action="store", dest='featureflag', default='CMAT',
							choices=['CMAT', 'aSLATM', 'FCHL', 'ACSF', 'BCAI'])
	parser.add_argument('--cutoff', help='cutoff distance for features',
							action="store", dest='cutoff', default=5.0)

	parser.add_argument('--targetflag', help='parameter to test',
							action="store", dest='targetflag', default='CCS')

	# Optional argments for train
	parser.add_argument('--searchmethod', help='Method for hyper-parameter search',
							action="store", dest='searchmethod', default='grid',
							choices=['grid', 'gaussian', 'random'])
	parser.add_argument('--param_ranges', help='Dictionary of parameter ranges for HPS',
							action="store", dest='param_ranges', default={})
	parser.add_argument('--param_logs', help='Dictionary (truth values only) for each parameter specifying whether parameter is on a log scale',
							action="store", dest='param_logs', default={})
	parser.add_argument('--cv_steps', help='Number of cross validations to perform',
	 						action="store", dest='cv_steps', default=5)
	parser.add_argument('--epochs', help='Number of HPS iterations to perform',
							action="store", dest='epochs', default=200)
	parser.add_argument('--logfile', help='Name of output log file',
							action="store", dest='logfile', default='')


	# Optional system arguments
	parser.add_argument('--python_env', help='Name of python environment to be used',
							action="store", dest='python_env', default='env_IMP')
	parser.add_argument('--system', help='System currently running',
	 						action="store", dest='system', default='localbox')
	parser.add_argument('--memory', help='Memory needed in submission scripts',
	 						action="store", dest='memory', default=3)
	parser.add_argument('--processors', help='Processors needed in submission scripts',
	 						action="store", dest='processors', default=1)
	parser.add_argument('--walltime', help='Walltime required for submission',
	 						action="store", dest='walltime', default='100:00:00')

	# Optional argments for predict


	# Optional argments for grid search
	parser.add_argument('--grid_density', help='Point density for grid search',
	 						action="store", dest='grid_density', default='localbox')

	# Optional arguments for gaussian search
	parser.add_argument('--kappa', help='Kappa value for gaussian process HPS',
	 						action="store", dest='kappa', default=5.0)
	parser.add_argument('--xi', help='Xi value for gaussian process HPS',
							action="store", dest='xi', default=0.1)
	parser.add_argument('--random', help='Frequency of random selection in gaussian process search',
							action="store", dest='random', default=0)


	# Parse arguments into args object
	args = parser.parse_args()

	# check param flag:
	param = targetflag.flag_to_target(args.targetflag)
	if param == 0:
		print('Invalid parameter flag')
		sys.exit(0)
	# check flag combination
	if not flag_combos.check_combination(args.modelflag, args.featureflag):
		print('Invalid model and feature combination: ', args.modelflag, args.featureflag)
		sys.exit(0)

	# Print pretty banner
	print_header_IMP()

	# set up submission file for model training
	if args.Command == "setup_train":

		CMD_trainmodel.setup_trainmodel(args)


	if args.Command == "train":
		if args.logfile == '':
			args.logfile = args.modelflag + '_' + args.featureflag + '_' + args.targetflag + '_' + args.searchmethod + '.log'
		CMD_trainmodel.trainmodel(args)


	elif args.Command == "predict":
		print('Not done yet')



	elif args.Command == 'compare':
		print('Not done yet')


	elif args.Command == 'test':
		print('Not done yet')































##

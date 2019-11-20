# master script to do everything IMPRESSION-y

import argparse
import sys
import json
sys.path.append('aE_lib/')

from util.flag_handler import hdl_targetflag, flag_combos
from util.argparse_wizard import run_wizard

from top_level import CMD_trainmodel, CMD_predict, CMD_compare

from util.header import print_header_IMP

import ast

import sys

def trace(frame, event, arg):
	print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
	return trace

#sys.settrace(trace)




if __name__ == "__main__":
	# Argparser
	parser = argparse.ArgumentParser(description='IMPRESSION')
	# Define command to perform
	parser.add_argument('Command', help='IMPRESSION command',
							choices=['setup_train', 'train', 'setup_predict', 'predict', 'compare', 'test'])

	parser.add_argument('--prefs', help='How to obtain settings for command, leave blank to use command line input, specify json file, or <wizard> to run wizard',
							 default='')

	parser.add_argument('--training_set', help='Training dataset file(s), either single csv/pkl file or one/multiple nmredata files',
						 action="store", dest='training_set', default='None')


	parser.add_argument('--store_datasets', help='Option to store datasets as pickle files for later use',
						 action="store", dest='store_datasets', default='False')


	parser.add_argument('--target_list', help='Optional list of targets to go through',
						 action="store", dest='target_list', default=[])
	parser.add_argument('--targetflag', help='parameter to test',
							action="store", dest='targetflag', default='CCS')

	# Optional arguments for train/test
	parser.add_argument('--modelflag', help='type of model to use',
							action="store", dest='modelflag', default='KRR',
							choices=['KRR', 'FCHL', 'NN', 'TFM'])
	parser.add_argument('--featureflag', help='type of features to use',
							action="store", dest='featureflag', default='CMAT',
							choices=['CMAT', 'aSLATM', 'FCHL', 'ACSF', 'BCAI'])
	parser.add_argument('--cutoff', help='cutoff distance for features',
							action="store", dest='cutoff', default=5.0)

	# Optional argments for train
	parser.add_argument('--searchmethod', help='Method for hyper-parameter search',
							action="store", dest='searchmethod', default='grid',
							choices=['grid', 'gaussian', 'random'])
	parser.add_argument('--feature_optimisation', help='HPS includes feature parameters',
							action="store", dest='feature_optimisation', default='True',
							choices=['True', 'False'])
	parser.add_argument('--feature_file', help='File containing pre-made feature dataset object',
							action="store", dest='feature_file', default='none')

	parser.add_argument('--param_ranges', help='Dictionary of parameter ranges for HPS',
							action="store", dest='param_ranges', default='{}')
	parser.add_argument('--param_logs', help='Dictionary (truth values only) for each parameter specifying whether parameter is on a log scale',
							action="store", dest='param_logs', default='{}')
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
	parser.add_argument('--var', help='Number of pre-prediction variance models',
						action="store", dest='var', default=0)
	parser.add_argument('--models', help='Existing model(s) to use, list',
						 action="store", dest='models', default=['None'])
	parser.add_argument('--test_sets', help='Testing dataset(s), either lists of file search patterns or individual files, list',
						 action="store", dest='test_sets', default='None')


	# Optional argments for grid search
	parser.add_argument('--grid_density', help='Point density for grid search',
							 action="store", dest='grid_density', default=10)

	# Optional arguments for gaussian search
	parser.add_argument('--kappa', help='Kappa value for gaussian process HPS',
							 action="store", dest='kappa', default=5.0)
	parser.add_argument('--xi', help='Xi value for gaussian process HPS',
							action="store", dest='xi', default=0.1)
	parser.add_argument('--random', help='Frequency of random selection in gaussian process search',
							action="store", dest='random', default=0)

	# Optional arguments for make features
	parser.add_argument('--input_files', help='Files to create features from',
							 action="store", dest='input_files', default='none')
	parser.add_argument('--output_files', help='Files to store features in',
							 action="store", dest='output_files', default='none')


	# Parse arguments into args object
	args = vars(parser.parse_args())

	COMMAND = args['Command']
	if args['prefs'] == 'wizard':
		args = run_wizard(args)
		pref_file = 'IMPRESSION_settings.json'
		args['prefs'] = pref_file
		json.dump(args, open(pref_file, 'w'), indent=4)

		with open(pref_file, 'r') as fp:
			args = json.load(fp)

	elif len(args['prefs']) == 0:
		pass
	else:
		print('Reading settings from file ', args['prefs'])
		try:
			with open(args['prefs'], 'r') as json_file:
				args = json.load(json_file)
		except Exception as E:
			print('Error reading preferences file ', file)
			print(E)
			print('Exiting. . .')
			sys.exit(0)

	args['Command'] = COMMAND
	# check param flag:
	param = hdl_targetflag.flag_to_target(args['targetflag'])
	if param == 0:
		print('Invalid parameter flag')
		sys.exit(0)
	# check flag combination
	if not flag_combos.check_combination(args['modelflag'], args['featureflag']):
		print('Invalid model and feature combination: ', args['modelflag'], args['featureflag'])
		sys.exit(0)

	# Print pretty banner
	print_header_IMP()

	# set up submission file for model training
	if args['Command'] == "setup_train":

		if len(target_list) > 0:
			for target in target_list:
				args['targetflag'] = target
				CMD_trainmodel.setup_trainmodel(args)
		else:
			CMD_trainmodel.setup_trainmodel(args)

		print('Training submission file created . . .')

	if args['Command'] == "train":
		if args['logfile'] == '':
			args['logfile'] = args['modelflag'] + '_' + args['featureflag'] + '_' + args['targetflag'] + '_' + args['searchmethod'] + '.log'
		CMD_trainmodel.trainmodel(args)


	elif args['Command'] == "setup_predict":
		CMD_predict.setup_predict(args)

	elif args['Command'] == "predict":
		CMD_predict.predict(args)



	elif args['Command'] == 'compare':
		print('Not done yet')


	elif args['Command'] == 'test':
		print('Not done yet')


























##

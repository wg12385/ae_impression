# Copyright 2020 Will Gerrard
#This file is part of autoENRICH.

#autoENRICH is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoENRICH is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoENRICH.  If not, see <https://www.gnu.org/licenses/>.

# Argument parser
import argparse
# To quit things
import sys
# Preferences files are all in json
import json

# Functions for checking input flags
from autoENRICH.util.flag_handler import hdl_targetflag, flag_combos
# Preferences wizard function
from autoENRICH.util.argparse_wizard import run_wizard, get_minimal_args
# Import main command functions
from autoENRICH.top_level import CMD_trainmodel, CMD_predict
# Import pretty banner printing function (for ego purposes only)
from autoENRICH.util.header import print_header_IMP
# Used for memory and code tracing
import ast
import tracemalloc
import cProfile
import pstats
from pstats import SortKey

# Define tracer, used to trace code execution line by line
def trace(frame, event, arg):
	print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
	return trace

if __name__ == "__main__":
	# Argparser
	parser = argparse.ArgumentParser(description='IMPRESSION')
	# Define command to perform
	parser.add_argument('Command', help='IMPRESSION command',
							choices=['setup_train', 'train', 'setup_predict', 'predict', 'test'])
	# Get preferences
	parser.add_argument('--prefs', help='How to obtain settings for command, leave blank to use command line input, specify json file, or <wizard> to run wizard',
							 default='')
	# Define training set, only needed for train and setup_train commands
	parser.add_argument('--training_set', help='Training dataset file(s), either single csv/pkl file or one/multiple nmredata files',
						 action="store", dest='training_set', default='')
	# Define output directory for non log files
	parser.add_argument('--output_dir', help='Output directory for non log files',
						action="store", dest='output_dir', default='')
	# Check whether to store prepared datasets or not
	parser.add_argument('--store_datasets', help='Option to store datasets as pickle files for later use',
						 action="store", dest='store_datasets', default='')

	# Define list of target parameters
	parser.add_argument('--target_list', help='Optional list of targets to go through',
						 action="store", dest='target_list', default=[])
	# Define single target
	parser.add_argument('--targetflag', help='target NMR parameter',
							action="store", dest='targetflag', default='')

	# Optional arguments for train/test
	# Type of model to use
	parser.add_argument('--modelflag', help='type of model to use',
							action="store", dest='modelflag', default='',
							choices=['KRR', 'FCHL', 'NN', 'TFM'])
	# Type of features to user
	parser.add_argument('--featureflag', help='type of features to use',
							action="store", dest='featureflag', default='',
							choices=['CMAT', 'aSLATM', 'FCHL', 'ACSF', 'BCAI', ''])

	# Optional argments for train
	# Method for hyper-parameter search
	parser.add_argument('--searchflag', help='Method for hyper-parameter search',
							action="store", dest='searchflag', default='',
							choices=['grid', 'gaussian', 'random', ''])
	# Whether to optimise features in hyper-parameter optimisation (uses defaults for features otherwise)
	# This offers a pretty massive speed improvement especially for some features (aSLATM)
	parser.add_argument('--feature_optimisation', help='HPS includes feature parameters',
							action="store", dest='feature_optimisation', default='',
							choices=['True', 'False', ''])
	# Can supply pre-made features to use
	parser.add_argument('--feature_file', help='File containing pre-made feature dataset object',
							action="store", dest='feature_file', default='')
	# Define ranges for model parameters to optimise
	parser.add_argument('--param_ranges', help='Dictionary of parameter ranges for HPS',
							action="store", dest='param_ranges', default='{}')
	# Specify whether to treat paramater is linear of log scale (helps with optimisation)
	parser.add_argument('--param_logs', help='Dictionary (truth values only) for each parameter specifying whether parameter is on a log scale',
							action="store", dest='param_logs', default='{}')
	# Number of cross validation steps to perform
	parser.add_argument('--cv_steps', help='Number of cross validations to perform',
							 action="store", dest='cv_steps', default=-404)
	# How many HPS iterations to perform
	parser.add_argument('--epochs', help='Number of HPS iterations to perform',
							action="store", dest='epochs', default=-404)
	# Specify output logfile
	parser.add_argument('--logfile', help='Name of output log file',
							action="store", dest='logfile', default='')


	# Optional system arguments
	parser.add_argument('--python_env', help='Name of python environment to be used',
							action="store", dest='python_env', default='')
	# Define system (allows for HPC cluster specific submission)
	parser.add_argument('--system', help='System currently running',
							 action="store", dest='system', default='')
	# How much memory to allow in a HPC submission
	parser.add_argument('--memory', help='Memory needed in submission scripts',
							 action="store", dest='memory', default=-404)
	# How many processors to request in HPC submission
	parser.add_argument('--processors', help='Processors needed in submission scripts',
							 action="store", dest='processors', default=-404)
	# How much walltime to request in HPC submission
	parser.add_argument('--walltime', help='Walltime required for submission',
							 action="store", dest='walltime', default='')

	# Optional argments for predict
	# How many variance models to use
	parser.add_argument('--var', help='Number of pre-prediction variance models',
						action="store", dest='var', default=-404)
	# Define existing models to use for predictions
	parser.add_argument('--models', help='Existing model(s) to use, list',
						 action="store", dest='models', default=[])
	# Current prediction model, not sure if this is used
	parser.add_argument('--model', help='Current prediction model',
						action="store", dest='model', default='none')
	# Define set(s) of molecules to perform predictions on
	parser.add_argument('--test_sets', help='Testing dataset(s), either lists of file search patterns or individual files, list',
						 action="store", dest='test_sets', default=[])


	# Optional argments for grid search
	# Define density of points in grid
	parser.add_argument('--grid_density', help='Point density for grid search',
							 action="store", dest='grid_density', default=-404)

	# Optional arguments for gaussian search
	parser.add_argument('--kappa', help='Kappa value for gaussian process HPS',
							 action="store", dest='kappa', default=-404)
	parser.add_argument('--xi', help='Xi value for gaussian process HPS',
							action="store", dest='xi', default=-404)
	parser.add_argument('--random', help='Frequency of random selection in gaussian process search',
							action="store", dest='random', default=-404)
	parser.add_argument('--load', help='Load previous search from logfile before starting',
							action="store", dest='load', default=False)

	# Optional arguments for make features
	parser.add_argument('--input_files', help='Files to create features from',
							 action="store", dest='input_files', default='none')
	parser.add_argument('--output_files', help='Files to store features in',
							 action="store", dest='output_files', default='none')

	# Be very careful with this, having a different size between the training and
	# test sets screws up the predictions
	parser.add_argument('--max_size', help='Maximum molecule size',
							action="store", dest='max_size', default=-404)
	# Option to trace the code execution, for bug hunting
	parser.add_argument('--tracecode', help='Trace the code execution',
							action="store_true", dest='tracecode', default=False)
	# Option to trace memory usage, prints the biggest 10 memory blocks at the end of execution
	parser.add_argument('--tracemem', help='Trace the memory usage',
								action="store_true", dest='tracemem', default=False)
	# Option to trace code execution by time
	parser.add_argument('--tracetime', help='Trace the execution timings',
								action="store_true", dest='tracetime', default=False)


	# Parse arguments into args object
	args = vars(parser.parse_args())

	# Optional, trace code execution
	if args['tracecode']:
		sys.settrace(trace)
	# Optional, trace memory
	if args['tracemem']:
		tracemalloc.start()

	if args['tracetime']:
		TRACETIME = True
		pr = cProfile.Profile()
		pr.enable()
	else:
		TRACETIME = False

	user_args = args
	default_args = get_minimal_args()
	file_args = {}
	wiz_args = {}

	# Preserve command argument whilst messing about with preferences / args
	COMMAND = args['Command']
	# Run preferences wizard to ask user for preference choices
	if args['prefs'] == 'wizard':
		wiz_args = run_wizard(args)
		pref_file = 'settings_' + args['modelflag'] + '_' + args['featureflag'] + '_' + args['targetflag'] + '.json'
		args['prefs'] = pref_file
		json.dump(args, open(pref_file, 'w'), indent=4)
		with open(pref_file, 'r') as fp:
			args = json.load(fp)
	# Run preferences wizard but select all default options (minimal user input)
	elif args['prefs'] == 'default':
		wiz_args = run_wizard(args, default=True)
		pref_file = 'settings_' + args['modelflag'] + '_' + args['featureflag'] + '_' + args['targetflag'] + '.json'
		args['prefs'] = pref_file
		json.dump(args, open(pref_file, 'w'), indent=4)
		with open(pref_file, 'r') as fp:
			args = json.load(fp)
	# Use default preferences for every argument and generate a settings file to edit
	elif args['prefs'] == 'generate':
		pref_file = 'IMPRESSION_settings.json'
		args['prefs'] = pref_file
		json.dump(args, open(pref_file, 'w'), indent=4)
		print('Template preferences file generated')
		sys.exit()
	# Else read preferences file or shout at the user for not specifiying one
	else:
		print('Reading settings from file ', args['prefs'])
		try:
			with open(args['prefs'], 'r') as json_file:
				file_args = json.load(json_file)
		except Exception as E:
			print('Error reading preferences file ', args['prefs'])
			print('You must specify a preferences file, or generate one using --prefs default')
			print(E)
			print('Exiting. . .')
			sys.exit(0)

	for arg in user_args:
		user = False
		if type(user_args[arg]) is int:
			if user_args[arg] != -404:
				user = True
		elif type(user_args[arg]) is str:
			if user_args[arg] not in ['', 'none', '{}']:
				user = True
		elif type(user_args[arg]) is dict:
			if len(user_args[arg]) != 0:
				user = True
		elif type(user_args[arg]) is list:
			if len(user_args[arg]) != 0:
				user = True

		file = False
		if arg in file_args:
			if type(file_args[arg]) is int:
				if file_args[arg] != -404:
					file = True
			elif type(file_args[arg]) is str:
				if file_args[arg] not in ['', 'none', '{}']:
					file = True
			elif type(file_args[arg]) is dict:
				if len(file_args[arg]) != 0:
					file = True
			elif type(file_args[arg]) is list:
				if len(file_args[arg]) != 0:
					file = True


		if user:
			args[arg] = user_args[arg]
		elif file:
			args[arg] = file_args[arg]
		elif arg in wiz_args:
			args[arg] = wiz_args[arg]
		elif arg in default_args:
			args[arg] = default_args[arg]
		else:
			continue

	# Restore command argument
	args['Command'] = COMMAND

	# Unless making predictions check combination of feature / model
	if args['Command'] not in ['predict', 'setup_predict']:
		# check target flag is valid (nJxy or XCS):
		target = hdl_targetflag.flag_to_target(args['targetflag'])
		# 0 is the bad number
		if target == 0:
			print('Invalid target flag, ', args['targetflag'])
			sys.exit(0)

		# check flag combination for feature / model
		if not flag_combos.check_combination(args['modelflag'], args['featureflag']):
			print('Invalid model and feature combination: ', args['modelflag'], args['featureflag'])
			sys.exit(0)

	# Print pretty banner
	print_header_IMP()

	# set up submission file for model training
	if args['Command'] == "setup_train":
		# If multiple targets, setup submission script for each target
		if len(args['target_list']) > 0:
			# Loop through targets
			for target in args['target_list']:
				# Assign targetflag
				args['targetflag'] = target
				# Define preferences file
				pref_file = 'settings_HPS_' + args['modelflag'] + '_' + args['featureflag'] + '_' + args['targetflag'] + '_' + args['searchflag'] + '.json'
				# Assign to args dict
				args['prefs'] = pref_file
				# Dump prefs in pref file
				json.dump(args, open(pref_file, 'w'), indent=4)
				# Run train_model setup function to create submission file
				CMD_trainmodel.setup_trainmodel(args)
		else:
			# Run train_model setup function to create submission file
			CMD_trainmodel.setup_trainmodel(args)
		# yay success. . .
		print('Training submission file created . . .')

	# Train a model
	if args['Command'] == "train":
		# Check for defined logfile
		if args['logfile'] == '':
			# If no logfile defined, make a sensible one
			args['logfile'] = args['modelflag'] + '_' + args['featureflag'] + '_' + args['targetflag'] + '_' + args['searchflag'] + '.log'

		# Tell the user what model they are training
		print('Training model: ', args['modelflag'] , args['featureflag'], args['targetflag'], args['searchflag'])
		# Train the model
		CMD_trainmodel.trainmodel(args)

	# Set up a prediction
	# This seems like a dumb thing to have, but it facilitates predictions on HPC clusters,
	# and its easier to edit a preferences file than it is to define all the flags
	elif args['Command'] == "setup_predict":
		# Set these flags to empty because these are taken from the specified model objects
		args['targetflag'] = ''
		args['modelflag'] = ''
		args['featureflag'] = ''
		# Define basic preferences file if none specified
		if args['prefs'] == '':
			pref_file = 'settings_predict.json'
			args['prefs'] = pref_file
		json.dump(args, open(args['prefs'], 'w'), indent=4)
		CMD_predict.setup_predict(args)

	# Make predictions from a molecule
	elif args['Command'] == "predict":
		CMD_predict.predict(args)

	# Do code testing, or in this case print sarcastic message
	elif args['Command'] == 'test':
		print('Not done yet, if will wasnt so lazy we would have some nice test code ')

	# Output for memory trace
	try:
		if args['tracemem']:
			snapshot = tracemalloc.take_snapshot()
			top_stats = snapshot.statistics('lineno')

			print("[ Top 10 ]")
			for stat in top_stats[:10]:
			    print(stat)
	except:
		print('No memory trace option')

	#if args['tracetime']:
	if TRACETIME:
		pr.disable()
		ps = pstats.Stats(pr).sort_stats('time')
		ps.print_stats(10)

















##

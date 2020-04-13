

# Test setup_train and train commands
'''
from autoenrich.top_level import CMD_trainmodel
from autoenrich.util.argparse_wizard import get_minimal_args, run_wizard
from autoenrich.util.flag_handler import hdl_targetflag, flag_combos
from autoenrich.util.arguments.argparser import combine_args

import glob
import os
import sys
import traceback

def construct_testcases():

	argsets = []

	for targetflag in ['HCS', 'CCS', 'NCS', '1JCH', '3JHH']:
		for modelflag in ['KRR', 'FCHL']: # NN TFM
			for featureflag in ['CMAT', 'aSLATM', 'FCHL', 'ACSF', 'BCAI']:
				for searchflag in ['grid', 'gaussian', 'random']:
					args = get_minimal_args()
					args['epochs'] = 1
					args['cv_steps'] = 2
					args['grid_density'] = 1
					args['load'] = False
					args['targetflag'] = targetflag
					args['modelflag'] = modelflag
					args['featureflag'] = featureflag
					args['searchflag'] = searchflag
					args['output_dir'] = './test_models/'
					argsets.append(args)
	return argsets

def test_setup():
	status = 'Pass'
	cases = construct_testcases()
	for case in cases:
		case['prefs'] = 'default'
		case['Command'] = 'setup_train'
		case['training_set'] = 'test_files/train*.sdf'
		case['python_env'] = 'IMPgen1'

		CMD_trainmodel.setup_trainmodel(case)

	print('Test status: ', status)

	return status


def test_train():
	status = 'Pass'
	cases = construct_testcases()
	for case in cases:
		print('CASE _____', case['featureflag'], case['modelflag'], case['targetflag'], case['searchflag'])
		case['prefs'] = 'default'
		case['Command'] = 'train'
		case['training_set'] = 'test_files/train*.sdf'
		case['python_env'] = 'IMPgen1'
		case['feature_optimisation'] = 'True'

		wiz_args = run_wizard(case, default=True)
		case = combine_args(case, {}, wiz_args, case)
		print(case['epochs'])
		# check target flag is valid (nJxy or XCS):
		target = hdl_targetflag.flag_to_target(case['targetflag'])
		# 0 is the bad number
		if target == 0:
			print('Invalid target flag, ', case['targetflag'])
			continue

		# check flag combination for feature / model
		if not flag_combos.check_combination(case['modelflag'], case['featureflag']):
			print('Invalid model and feature combination: ', case['modelflag'], case['featureflag'])
			continue

		CMD_trainmodel.trainmodel(case)

	print('Test status: ', status)

	return status



if __name__ == "__main__":
	status1 = test_setup()
	submit_files = glob.glob('*.submit')
	for file in submit_files:
		os.remove(file)

	status2 = test_train()
	files_to_remove = glob.glob('*.pkl')
	files_to_remove.extend(glob.glob('*.log'))
	for file in files_to_remove:
		os.remove(file)
	print(status1, status2)
'''

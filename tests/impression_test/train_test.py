

# Test setup_train and train commands

from autoENRICH.top_level import CMD_trainmodel
from autoENRICH.util.argparse_wizard import get_minimal_args


def construct_testcases():

	argsets = []

	for targetflag in ['HCS', 'CCS', 'NCS', '1JCH', '3JHH']:
		for modelflag in ['KRR', 'FCHL', 'NN', 'TFM']:
			for featureflag in ['CMAT', 'aSLATM', 'FCHL', 'ACSF', 'BCAI']:
				for searchflag in ['grid', 'gaussian', 'random']:
					args = get_minimal_args()
					args['epochs'] = 1
					args['cv_steps'] = 2
					args['targetflag'] = targetflag
					args['modelflag'] = modelflag
					args['featureflag'] = featureflag
					args['searchflag'] = searchflag
					argsets.append(args)
	return argsets

def test_setup():
	status = 'Pass'
	cases = construct_testcases()
	for testcase in cases:

		case['Command'] = 'setup_train'
		case['training_set'] = 'test_files/train*.sdf'
		case['python_env'] = 'IMPgen1'

		try:
			CMD_trainmodel.setup_trainmodel(case)
		except Exception as e:
			print('ERROR IN TEST -------------------------------------------')
			print(case)
			print(e)
			status = 'Fail'

	print('Test status: ', status)

	return status


def test_train():
	status = 'Pass'
	cases = construct_testcases()
	for testcase in cases:

		case['Command'] = 'train'
		case['training_set'] = 'test_files/train*.sdf'
		case['python_env'] = 'IMPgen1'

		try:
			CMD_trainmodel.trainmodel(case)
		except Exception as e:
			print('ERROR IN TEST -------------------------------------------')
			print(case)
			print(e)
			status = 'Fail'

	print('Test status: ', status)

	return status



if __name__ == "__main__":
	status1 = test_setup()
	status2 = test_train()
	print(status1, status2)


# Tests setup_predict and predict functions
from autoENRICH.top_level import CMD_predict
from autoENRICH.util.argparse_wizard import get_minimal_args
import glob

def get_test_cases():

	test_cases = {}

		for modelflag in ['KRR', 'FCHL', 'NN', 'TFM']:
			for featureflag in ['CMAT', 'aSLATM', 'FCHL', 'ACSF', 'BCAI']:
				for searchflag in ['grid', 'gaussian', 'random']:
					args = get_minimal_args()
					args['var'] = 2
					for targetflag in ['HCS', 'CCS', 'NCS', '1JCH', '3JHH']:
						model = modelflag + '_' + targetflag + '_' + featureflag + '_' + searchflag + '.pkl'
						args['models'].append(model)
					test_cases.append(model)

	return test_cases

def test_setup():
	status = 'Pass'
	case = get_minimal_args
	case['Command'] = 'setup_predict'
	case['test_sets'] = ['test_files/test*.sdf']
	try:
		CMD_predict.setup_predict(case)
	except Exception as e:
		print('ERROR IN TEST -------------------------------------------')
		print(case)
		print(e)
		status = 'Fail'

	return status

def test_predict():
	status = 'Pass'

	cases = get_test_cases()

	for case in cases:
		case['Command'] = 'predict'
		case['test_sets'] = ['test_files/test*.sdf']
		try:
			CMD_predict.predict(case)
		except Exception as e:
			print('ERROR IN TEST -------------------------------------------')
			print(case)
			print(e)
			status = 'Fail'

	return status




if __name__ == "__main__":
	status1 = test_setup()
	status2 = test_train()
	print(status1, status2)

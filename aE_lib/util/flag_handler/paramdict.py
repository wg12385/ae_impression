


def construct_param_dict(modelflag, featureflag, targetflag):

	params = {}
	param_ranges = {}
	param_logs = {}

	if modelflag == 'KRR':
		param_ranges['sigma'] = [-5, 5]
		param_logs['sigma'] = 'log'

		param_ranges['lamda'] = [-10, 1]
		param_logs['lamda'] = 'log'

	elif modelflag == 'FCHL':
		param_ranges['sigma'] = [-5, 5]
		param_logs['sigma'] = 'log'

		param_ranges['lamda'] = [-10, 1]
		param_logs['lamda'] = 'log'

	elif modelflag == 'NN':
		# add some params
		pass

	elif modelflag == 'TFM':
		# add some params
		pass


	if featureflag == 'CMAT':
		param_ranges['cutoff'] = [0, 10]
		param_logs['cutoff'] = 'lin'

	elif featureflag == 'aSLATM':
		param_ranges['cutoff'] = [0, 10]
		param_logs['cutoff'] = 'lin'

	elif featureflag == 'FCHL':
		param_ranges['cutoff'] = [0, 10]
		param_logs['cutoff'] = 'lin'

	elif featureflag == 'ACSF':
		param_ranges['cutoff'] = [0, 10]
		param_logs['cutoff'] = 'lin'

	elif featureflag == 'BCAI':
		param_ranges['cutoff'] = [0, 10]
		param_logs['cutoff'] = 'lin'


	return param_ranges, param_logs






























###



def check_combination(modelflag, featureflag):


	if modelflag == 'KRR':
		if featureflag not in ['CMAT', 'aSLATM', 'ACSF']:
			return False

	elif modelflag == 'FCHL':
		if featureflag != 'FCHL':
			return False


	elif modelflag == 'NN':
		return True


	elif modelflag == 'TFM':
		if featureflag != 'BCAI':
			return False

	else:
		return False


	return True

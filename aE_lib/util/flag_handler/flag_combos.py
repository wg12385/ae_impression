

def check_combination(modelflag, featureflag):


	if modelflag == 'KRR':
		if featureflag != 'CMAT' and featureflag != 'aSLATM':
			return False

	elif modelflag == 'FCHL':
		if featureflag != 'FCHL':
			return False


	elif modelflag == 'NN':
		return True


	elif modelflag == 'TFM':
		if featureflag != 'TFM':
			return False

	else:
		return False


	return True

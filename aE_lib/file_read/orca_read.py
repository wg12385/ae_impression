

def get_opt_status(file):
	status = 'unknown'
	try:
		with open(file, 'r') as f:
			for line in f:
				if 'SUCCESS' in line:
					status = 'successful'
				if 'ERROR' in line:
					status = 'failed'
	except Exception as e:
		status = 'not submitted'

	return status

def get_nmr_status(file):
	status = 'unknown'
	with open(file, 'r') as f:
		for line in f:
			if 'SUCCESS' in line:
				status = 'successful'
			if 'ERROR' in line:
				status = 'failed'
	return status

def read_structure(file):



	return xyz, types, conn, coupling_len


def read_opt(file):
	# SCF Energy:    -1072.8219232141
	energy = 0.0
	with open(file ,'r') as f:
		for line in f:
			if 'SCF Energy:' in line:
				items=line.split()
				energy = float(items[-1])

	return energy

def read_nmr(file):



	return shift, coupling

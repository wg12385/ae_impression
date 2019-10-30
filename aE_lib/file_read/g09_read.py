
def get_opt_status(file):
	status = 'unknown'
	try:
		with open(filename, 'r') as f:
			for line in f:
				if 'SUCCESS' in line:
					status = 'successful'
				if 'ERROR' in line:
					status = 'failed'
	except:
		status = 'not submitted'

	return status

def get_nmr_status(file):
	status = 'unknown'
	with open(filename, 'r') as f:
		for line in f:
			if 'SUCCESS' in line:
				status = 'successful'
			if 'ERROR' in line:
				status = 'failed'
	return status

def read_structure(file):



	return xyz, types, conn, coupling_len


def read_opt(file):



	return self.xyz, energy

def read_nmr(file):



	return shift, coupling

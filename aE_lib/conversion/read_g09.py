import numpy as np


def chemread(filename, atomnumber):
	# Read in chemical shift information from gaussian log files
	# Values read in are actually chemical shielding tensors, need to be converted by processing function

	shift_array = np.zeros(atomnumber, dtype=np.float64)

	switch = False
	with open(filename, 'r') as f_handle:
		for line in f_handle:
			if "SCF GIAO Magnetic shielding tensor (ppm)" in line:
				switch = True
			if "Fermi Contact" in line:
				switch = False

			if switch:
				if "Isotropic" in line:
					items = line.split()
					try:
						num = int(items[0])
					except:
						continue

					shift_array[num-1] = float(items[4])
	# Return array is zero indexed
	return shift_array

def jread(filename, atomnumber):

	couplings = np.zeros((atomnumber, atomnumber), dtype=float)
	with open(filename) as f:
		skip = True
		for line in f:
			if "Total nuclear spin-spin coupling J (Hz):" in line:
				skip = False
				continue
			elif "End of Minotr" in line:
				skip = True
				continue
			elif skip:
				continue

			if "D" not in line:
				tokens = line.split()
				i_indices = np.asarray(tokens, dtype=int)
			else:
				tokens = line.split()
				index_j = int(tokens[0]) - 1
				for i in range(len(tokens)-1):
					index_i = i_indices[i] - 1
					coupling = float(tokens[i+1].replace("D","E"))
					couplings[index_i, index_j] = coupling
					couplings[index_j, index_i] = coupling
	return couplings

def eread(filename):
	energy = 404.404
	print('FUNCTION NOT WRITTEN YET, eread, read_g09.py')

	return energy

def check_opt_status(filename):
	status = 'unknown'
	with open(filename, 'r') as f:
		for line in f:
			if 'Normal termination' in line:
				status = 'successful'
			if 'Error termination' in line:
				status = 'failed'

	return status



def check_nmr_status(filename):
	status = 'unknown'
	with open(filename, 'r') as f:
		for line in f:
			if 'Normal termination' in line:
				status = 'successful'
			if 'Error termination' in line:
				status = 'failed'

	return status

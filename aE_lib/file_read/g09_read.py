import numpy as np


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


def read_functional(file):
	functional = ''
	basisset = ''
	with open(file, 'r') as f:
		for line in f:
			# assumes nmrcommand of format:
			# #T nmr(stuff)functional/basisset stuff...
			if '#T nmr' in line:
				try:
					items = line.split()
					nmrcommand = items[2]
					functional = nmrcommand.split('/')[0]
					basisset = nmrcommand.split('/')[1]
				
					break
				except Exception as e:
					print(line, items)
					print(e)

	return functional, basisset


def read_nmr(file, atomnumber):

	shift_array = np.zeros(atomnumber, dtype=np.float64)

	switch = False
	with open(file, 'r') as f_handle:
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


	couplings = np.zeros((atomnumber, atomnumber), dtype=float)
	with open(file, 'r') as f:
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
					couplings[index_i][index_j] = coupling
					couplings[index_j][index_i] = coupling

	return shift_array, couplings

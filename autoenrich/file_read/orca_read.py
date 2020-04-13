# Copyright 2020 Will Gerrard
#This file is part of autoenrich.

#autoenrich is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoenrich is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import os

# Get status of ORCA optimisation
def get_opt_status(file):
	# Input:
	#	file: filename

	# Returns: status

	if not os.path.isfile(file):
		return 'not_submitted'
	else:
		status = 'unknown'
		finished = False
		with open(file, 'r') as f:
			for line in f:
				if 'SUCCESS' in line:
					status = 'successful'
				if 'ERROR' in line:
					status = 'failed'
				if '****ORCA TERMINATED NORMALLY****' in line:
					finished = True
		if finished:
			return status
		else:
			return 'unknown'

# Get status of ORCA NMR calculation
def get_nmr_status(file):
		# Input:
		#	file: filename

		# Returns: status

	if not os.path.isfile(file):
		return 'not submitted'
	else:
		status = 'unknown'
		finished = False
		with open(file, 'r') as f:
			for line in f:
				if 'SUCCESS' in line:
					status = 'successful'
				if 'ERROR' in line:
					status = 'failed'
				if '****ORCA TERMINATED NORMALLY****' in line:
					finished = True
		if finished or status == 'failed':
			return status
		else:
			return 'unknown'

# Read Optimisation energy from ORCA optimisation file
def read_opt(file):
	# Input:
	#	file: filename
	# Returns: energy

	# Go through file looking for this line
	# SCF Energy:    -1072.8219232141
	energy = 0.0
	with open(file ,'r') as f:
		for line in f:
			if 'FINAL SINGLE POINT ENERGY' in line:
				items=line.split()
				energy = float(items[-1])

	return energy

def read_functional(file):

	functional = ''
	basisset = ''

	with open(file, 'r') as f:
		for line in f:
			if '|  1> !' in line:
				functional = line.split()[3]
				basisset = line.split()[4]

	return functional, basisset


# Read NMR information from ORCA NMR log files
def read_nmr(file):

	shiftswitch = False
	shifts = []
	cplswitch = False
	cpls = []
	strucswitch = False
	atoms = 0

	with open(file ,'r') as f:
		for line in f:

			if 'CARTESIAN COORDINATES (A.U.)' in line:
				strucswitch = True

			if 'CHEMICAL SHIELDING SUMMARY (ppm)' in line:
				shiftswitch = True
				cplswitch = False
				strucswitch = False

			if 'NMR SPIN-SPIN COUPLING CONSTANTS' in line:
				shiftswitch = False
				strucswitch = False
				cplswitch = True

			items = line.split()
			if len(items) == 0:
				continue

			if strucswitch and len(items) == 8 and items[0] != 'NO':
				try:
					atoms = int(items[0]) + 1
				except:
					continue

			if shiftswitch and len(items) == 4:
				try:
					int(items[0])
					float(items[2])
					float(items[3])
				except:
					continue

				shifts.append([int(items[0]), float(items[2])])

			if cplswitch and len(items) in [6, 10]:
				if 'NUCLEUS A' in line and len(items) == 10:
					a = int(items[4])
					b = int(items[9])

				if items[0] == 'Total' and len(items) == 6:
					c = float(items[5])

					cpls.append([a, b, c])

	shift = np.zeros((atoms), dtype=np.float64)
	for sh in shifts:
		shift[sh[0]] = sh[1]

	coupling = np.zeros((atoms, atoms), dtype=np.float64)
	for cp in cpls:
		coupling[cp[0]][cp[1]] = cp[2]
		coupling[cp[1]][cp[0]] = cp[2]

	return shift, coupling



















##

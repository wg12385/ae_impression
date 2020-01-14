# Copyright 2020 Will Gerrard
#This file is part of autoENRICH.

#autoENRICH is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoENRICH is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with autoENRICH.  If not, see <https://www.gnu.org/licenses/>.


import numpy as np

def read_nmr(file, atoms):

	shift_array = np.zeros(atoms, dtype=np.float64)
	shift_var = np.zeros(atoms, dtype=np.float64)
	coupling_array = np.zeros((atoms, atoms), dtype=np.float64)
	coupling_len = np.zeros((atoms, atoms), dtype=np.int64)
	coupling_var = np.zeros((atoms, atoms), dtype=np.float64)

	with open(file, 'r') as f:
		shift_switch = False
		cpl_switch = False
		for line in f:
			if '> <NMREDATA_ASSIGNMENT>' in line:
				shift_switch = True

			if '> <NMREDATA_J>' in line:
				shift_switch = False
				cpl_switch = True


			if shift_switch:
				#  0    , -33.56610000   , 8    , 0.00000000     \
				items = line.split()

				try:
					int(items[0])
				except:
					continue

				shift_array[int(items[0])] = float(items[2])
				shift_var[int(items[0])] = float(items[6])

			if cpl_switch:
				#  0         , 4         , -0.08615310    , 3JON      , 0.00000000
				# ['0', ',', '1', ',', '-0.26456900', ',', '5JON', ',', '0.00000000']
				items = line.split()
				try:
					int(items[0])
				except:
					continue
				length = int(items[6].strip()[0])
				coupling_array[int(items[0])][int(items[2])] = float(items[4])
				coupling_var[int(items[0])][int(items[2])] = float(items[8])
				coupling_len[int(items[0])][int(items[2])] = length
				coupling_len[int(items[2])][int(items[0])] = length

	return shift_array, shift_var, coupling_array, coupling_var, coupling_len

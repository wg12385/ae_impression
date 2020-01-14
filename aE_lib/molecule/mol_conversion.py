# Copyright 2020 Will Gerrard
#This file is part of autoENRICH.

#autoENRICH is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoENRICH is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoENRICH.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from .convert import nmrmol_to_xyz
from .pybel_plugins import *

def nmrmol_to_pybmol(nmrmol):

	file = nmrmol.molid + 'tmp.xyz'

	nmrmol_to_xyz(nmrmol, file)

	pybmol = next(pyb.readfile('xyz', file))

	os.remove(file)

	return pybmol, nmrmol.molid

def pybmol_to_nmrmol(mol, molid):

	type_list, types = mol_read_type(mol)
	xyz = mol_read_xyz(mol)
	dist = mol_read_dist(mol)
	shift = np.zeros(len(types), dtype=np.float64)

	j_array = np.zeros((len(types),len(types)), dtype=np.float64)

	b1_paths = []
	b2_paths = []
	b3_paths = []

	b1_coupling = []
	b2_coupling = []
	b3_coupling = []

	for x in range(len(mol.atoms)):
		for y in range(x+1, len(mol.atoms)):
			b1_paths.extend(mol_find_all_paths(mol, x, y, 1))
			b2_paths.extend(mol_find_all_paths(mol, x, y, 2))
			b3_paths.extend(mol_find_all_paths(mol, x, y, 3))

	for path in b1_paths:
		if j_array[int(path[0])][int(path[1])] == 0:
			j = j_array[int(path[1])][int(path[0])]
		else:
			j = j_array[int(path[0])][int(path[1])]
		if int(path[0]) > int(path[1]):
			continue
		b1_coupling.append([int(path[0]), int(path[1]), 1, j])

	for path in b2_paths:
		if j_array[int(path[0])][int(path[1])] == 0:
			j = j_array[int(path[1])][int(path[0])]
		else:
			j = j_array[int(path[0])][int(path[1])]
		if int(path[0]) > int(path[2]):
			continue
		b2_coupling.append([int(path[0]), int(path[1]), int(path[2]), 2, j])

	for path in b3_paths:
		if j_array[int(path[0])][int(path[1])] == 0:
			j = j_array[int(path[1])][int(path[0])]
		else:
			j = j_array[int(path[0])][int(path[1])]
		if int(path[0]) > int(path[3]):
			continue
		b3_coupling.append([int(path[0]), int(path[1]), int(path[2]), int(path[3]), 3, j])


	shift_var = np.zeros(len(shift), dtype=np.float64)
	var1b = np.zeros(len(b1_coupling), dtype=np.float64)
	var2b = np.zeros(len(b2_coupling), dtype=np.float64)
	var3b = np.zeros(len(b3_coupling), dtype=np.float64)

	return molid, types, xyz, dist, shift, shift_var, b1_coupling, var1b, b2_coupling, var2b, b3_coupling, var3b

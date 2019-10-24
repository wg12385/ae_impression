# gaussian conversion functions

import numpy as np
import pybel as pyb
from .pybel_plugins import *
from .read_g09 import *
from reference.periodic_table import Get_periodic_table
from reference.tantillo import Get_tantillo_factors

def optg09_to_nmrdata(file, molid=0):

	scaling_factors = Get_tantillo_factors()

	mol = next(pyb.readfile('g09', file))
	type_list, types = mol_read_type(mol)
	xyz = mol_read_xyz(mol)
	dist = mol_read_dist(mol)
	shift = chemread(file, len(types))
	for i in range(len(shift)):
		if types[i] >= len(scaling_factors):
			continue

		shift[i] = (shift[i] - scaling_factors[types[i]][1]) / float(scaling_factors[types[i]][0])

	j_array = jread(file, len(types))

	molid = file.split('/')[-1].split('_')[1]

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


	return molid, types, xyz, dist, shift, b1_coupling, b2_coupling, b3_coupling


def MINIMIZE_g09_to_nmrdata(file, molid=0):

	scaling_factors = Get_tantillo_factors()

	mol = next(pyb.readfile('g09', file))

	mol.localopt(forcefield='mmff94', steps=2000)

	type_list, types = mol_read_type(mol)
	xyz = mol_read_xyz(mol)
	dist = mol_read_dist(mol)
	shift = chemread(file, len(types))
	for i in range(len(shift)):
		if types[i] >= len(scaling_factors):
			continue
		shift[i] = (shift[i] - scaling_factors[types[i]][1]) / float(scaling_factors[types[i]][0])
	j_array = jread(file, len(types))

	molid = file.split('/')[-1].split('_')[0]

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

def TENSOR_g09_to_nmrdata(file, molid=0):

	scaling_factors = Get_tantillo_factors()

	mol = next(pyb.readfile('g09', file))
	type_list, types = mol_read_type(mol)
	xyz = mol_read_xyz(mol)
	dist = mol_read_dist(mol)
	shift = chemread(file, len(types))

	j_array = jread(file, len(types))

	molid = file.split('/')[-1].split('_')[0]

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



def g09_to_nmrdata(file, molid=0):

	scaling_factors = Get_tantillo_factors()

	mol = next(pyb.readfile('g09', file))
	type_list, types = mol_read_type(mol)
	xyz = mol_read_xyz(mol)
	dist = mol_read_dist(mol)
	shift = chemread(file, len(types))
	for i in range(len(shift)):
		if types[i] >= len(scaling_factors):
			continue
		shift[i] = (shift[i] - scaling_factors[types[i]][1]) / float(scaling_factors[types[i]][0])
	j_array = jread(file, len(types))

	molid = file.split('/')[-1].split('_')[0] + file.split('/')[-1].split('_')[1]

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

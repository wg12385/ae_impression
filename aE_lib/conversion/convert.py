import numpy as np
import pybel as pyb
from .pybel_plugins import *
from .read_g09 import *
from reference.periodic_table import Get_periodic_table
from reference.tantillo import Get_tantillo_factors

from .cnv_gaussian import *
from .cnv_nxyz import *
from .cnv_syngenta import *
from .cnv_nmredata import *


import os


from .cnv_nmredata import nmrmol_to_nmredata

def labelmaker(cpl, mol):
	Periodic_table = Get_periodic_table()
	lent = len(cpl)-3
	label = str(lent) + str('J')
	if mol.types[int(cpl[0])] >= mol.types[int(cpl[lent])]:
		label = label + str(Periodic_table[mol.types[int(cpl[0])]]) + str(Periodic_table[mol.types[int(cpl[lent])]])
	else:
		label = label + str(Periodic_table[mol.types[int(cpl[lent])]]) + str(Periodic_table[mol.types[int(cpl[0])]])

	return label

def nmrmol_to_xyz(mol, outname, num=-404):
	periodic_table = Get_periodic_table()
	with open(outname, 'w') as f:
		print(len(mol.types), file=f)
		if num == -404:
			print(mol.molid, file=f)
		else:
			string = "{0:<10d}\t{1:<20s}".format(num, mol.molid)
			print(string, file=f)

		for i in range(len(mol.types)):
			string = "{i:<10s}\t{x:<10.6f}\t{y:<10.6f}\t{z:<10.6f}".format(i=periodic_table[mol.types[i]],
																			x=mol.xyz[i][0],
																			y=mol.xyz[i][1],
																			z=mol.xyz[i][2])
			print(string, file=f)

def pdb_to_nmrdummy(file):
	print(file)
	mol = next(pyb.readfile('pdb', file))
	type_list, types = mol_read_type(mol)
	xyz = mol_read_xyz(mol)
	dist = mol_read_dist(mol)
	shift = np.zeros(len(types), dtype=np.float64)

	j_array = np.zeros((len(types),len(types)), dtype=np.float64)

	#molid = file.split('/')[-1].split('_')[0]
	molid = file.split('/')[-1].split('.')[0]

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

def xyz_to_nmrdummy(file):

	mol = next(pyb.readfile('xyz', file))
	type_list, types = mol_read_type(mol)
	xyz = mol_read_xyz(mol)
	dist = mol_read_dist(mol)
	shift = np.zeros(len(types), dtype=np.float64)

	j_array = np.zeros((len(types),len(types)), dtype=np.float64)

	#molid = file.split('/')[-1].split('_')[0]
	molid = file.split('/')[-1].split('.')[0]

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

def mol2_to_nmrdummy(file):

	mol = next(pyb.readfile('mol2', file))
	type_list, types = mol_read_type(mol)
	xyz = mol_read_xyz(mol)
	dist = mol_read_dist(mol)
	shift = np.zeros(len(types), dtype=np.float64)

	j_array = np.zeros((len(types),len(types)), dtype=np.float64)

	#molid = file.split('/')[-1].split('_')[0]
	molid = file.split('/')[-1].split('.')[0]

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

def split_multiple_xyx(file):

	molid = file.split('/')[-1].split('.')[0]
	mol = next(pyb.readfile('xyz', file))
	type_list, types = mol_read_type(mol)
	xyz = mol_read_xyz(mol)

	structures = mol_splitmol(mol)

	periodic_table = Get_periodic_table()

	s = 0
	for structure in structures:
		if len(structure) < 5:
			continue
		s += 1
		with open(molid + str(s) + '.xyz', 'w') as f:
			print(len(structure), file=f)
			print(molid, file=f)
			for i in structure:
				string = "{i:<10s}\t{x:<10.6f}\t{y:<10.6f}\t{z:<10.6f}".format(i=periodic_table[types[i]],
																				x=xyz[i][0],
																				y=xyz[i][1],
																				z=xyz[i][2])
				print(string, file=f)

	print(molid, s)

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

def generic_pybel_read(file, type):
	mol = next(pyb.readfile(type, file))
	type_list, types = mol_read_type(mol)
	xyz = mol_read_xyz(mol)
	dist = mol_read_dist(mol)
	molid = file.split('/')[-1].split('.')[0]

	return molid, xyz, types, dist



###

import numpy as np
import pybel as pyb
from .pybel_plugins import *
from .read_g09 import *
from reference.periodic_table import Get_periodic_table
from reference.tantillo import Get_tantillo_factors

def labelmaker(cpl, mol):
	Periodic_table = Get_periodic_table()
	lent = len(cpl)-3
	label = str(lent) + str('J')
	if mol.types[int(cpl[0])] >= mol.types[int(cpl[lent])]:
		label = label + str(Periodic_table[mol.types[int(cpl[0])]]) + str(Periodic_table[mol.types[int(cpl[lent])]])
	else:
		label = label + str(Periodic_table[mol.types[int(cpl[lent])]]) + str(Periodic_table[mol.types[int(cpl[0])]])

	return label

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


def nxyz_to_nmrdata(name, path=''):
	with open(str(path)+str(name), 'r') as f:
		xyzswitch = False
		cplswitch = False
		l = 0

		for line in f:
			l += 1
			items = line.split()

			if l == 1:
				atoms = int(items[0])
				types = np.zeros(atoms, dtype=np.int64)
				xyz = np.zeros((atoms, 3), dtype=np.float64)
				dist = np.zeros((atoms, atoms), dtype=np.int64)
				shift = np.zeros(atoms, dtype=np.float64)
				continue

			if l == 2:
				try:
					molid = items[1]
				except:
					try:
						molid = items[0]
					except:
						molid = 0

			if 'END' in items:
				xyzswitch = False
				cplswitch = False

			if 'XYZ_DATA' in items:
				xyzswitch = True

			if 'COUPLING_DATA' in items:
				coupling1b = []
				coupling2b = []
				coupling3b = []
				cplswitch = True

			if xyzswitch:
				if len(items) != 6:
					continue
				i = int(items[0])
				types[i] = int(items[1])
				xyz[i][0] = float(items[2])
				xyz[i][1] = float(items[3])
				xyz[i][2] = float(items[4])
				shift[i] = float(items[5])

			if cplswitch:
				if len(items) < 4:
					continue
				try:
					int(items[1])
				except:
					continue

				if len(items) == 5:
					# (0)label at1 at2 len coupling
					coupling1b.append([int(items[1]), int(items[2]), int(items[3]), float(items[4])])
				if len(items) == 6:
					coupling2b.append([int(items[1]), int(items[2]), int(items[3]), int(items[4]), float(items[5])])
				if len(items) == 7:
					coupling3b.append([int(items[1]), int(items[2]), int(items[3]), int(items[4]), int(items[5]), float(items[6])])

	b1_coupling = np.asarray(coupling1b)
	b2_coupling = np.asarray(coupling2b)
	b3_coupling = np.asarray(coupling3b)

	return molid, types, xyz, dist, shift, b1_coupling, b2_coupling, b3_coupling

def nmrmol_to_nxyz(mol, outname, num=-404):
	with open(outname, 'w') as f:
		print(len(mol.types), file=f)
		if num == -404:
			print(mol.molid, file=f)
		else:
			string = "{0:<10d}\t{1:<20s}".format(num, mol.molid)
			print(string, file=f)

		print('', file=f)
		print('XYZ_DATA', file=f)
		for i in range(len(mol.types)):
			string = "{i:<10d}\t{typ:<10d}\t{x:<10.6f}\t{y:<10.6f}\t{z:<10.6f}\t{d:<10.6f}".format(i=i,
																						typ=int(mol.types[i]),
																						x=mol.xyz[i][0],
																						y=mol.xyz[i][1],
																						z=mol.xyz[i][2],
																						d=mol.shift[i])
			print(string, file=f)
		print('END', file=f)
		print('', file=f)
		print('COUPLING_DATA', file=f)

		for i in range(len(mol.coupling1b)):
			if mol.coupling1b[i][0] > mol.coupling1b[i][1]:
				continue
			label = labelmaker(mol.coupling1b[i], mol)
			string = "{lbl:<10s}\t{a1:<10d}\t{a2:<10d}\t{len:<10d}\t{coupling:<10.6f}".format(lbl=label,
																					a1=int(mol.coupling1b[i][0]),
																					a2=int(mol.coupling1b[i][1]),
																					len=int(mol.coupling1b[i][2]),
																					coupling=mol.coupling1b[i][3])
			print(string, file=f)

		for i in range(len(mol.coupling2b)):
			if mol.coupling2b[i][0] > mol.coupling2b[i][2]:
				continue
			label = labelmaker(mol.coupling2b[i], mol)
			string = "{lbl:<10s}\t{a1:<10d}\t{a2:<10d}\t{a3:<10d}\t{len:<10d}\t{coupling:<10.6f}".format(lbl=label,
																					a1=int(mol.coupling2b[i][0]),
																					a2=int(mol.coupling2b[i][1]),
																					a3=int(mol.coupling2b[i][2]),
																					len=int(mol.coupling2b[i][3]),
																					coupling=mol.coupling2b[i][4])
			print(string, file=f)

		for i in range(len(mol.coupling3b)):
			if mol.coupling3b[i][0] > mol.coupling3b[i][3]:
				continue
			label = labelmaker(mol.coupling3b[i], mol)
			string = "{lbl:<10s}\t{a1:<10d}\t{a2:<10d}\t{a3:<10d}\t{a4:<10d}\t{len:<10d}\t{coupling:<10.6f}".format(lbl=label,
																					a1=int(mol.coupling3b[i][0]),
																					a2=int(mol.coupling3b[i][1]),
																					a3=int(mol.coupling3b[i][2]),
																					a4=int(mol.coupling3b[i][3]),
																					len=int(mol.coupling3b[i][4]),
																					coupling=mol.coupling3b[i][5])
			print(string, file=f)
		print('END', file=f)

	return outname

def sygcml_to_nmrdata(name, path):

	mol = next(pyb.readfile('cml', str(path+name)))
	atoms = len(mol.atoms)
	xyz = mol_read_xyz(mol)
	type_list, types = mol_read_type(mol)
	shift = np.zeros(atoms, dtype=np.float64)
	dist = np.zeros((atoms,atoms), dtype=np.float64)

	b1_paths = []
	b2_paths = []
	b3_paths = []

	b1_coupling = []
	b2_coupling = []
	b3_coupling = []

	for x in set(types):
		for y in set(types):
			b1_paths.extend(mol_pathway_finder(mol, x, y, 1))
			b2_paths.extend(mol_pathway_finder(mol, x, y, 2))
			b3_paths.extend(mol_pathway_finder(mol, x, y, 3))

	j_array = np.zeros((atoms, atoms), dtype=np.float64)

	jch_array = np.zeros(atoms, dtype=np.float64)

	with open(str(path+name), 'r') as f:
		for line in f:
			items = line.split()
			if '<molecule id=' in line:
				molid = int(line.split('<molecule id="')[1].split('_')[0])
			if '<atom id=' in line:
				atomname = line.split('id="')[1].split('"')[0]
				type = line.split('elementType="')[1].split('"')[0]
				x = float(line.split('x3="')[1].split('"')[0])
				y = float(line.split('y3="')[1].split('"')[0])
				z = float(line.split('z3="')[1].split('"')[0])
			if '<scalar title="Exp1JCH"' in line:
			#if '<scalar title="NWCHEM_1JCH"' in line:
				try:
					coupling = float(line.split('dataType="xsd:string">')[1].split('</scalar>')[0])
				except Exception as e:
					print(name, e, line)
					continue

				for i in range(atoms):
					if x == xyz[i][0] and y == xyz[i][1] and z == xyz[i][2]:
						jch_array[i] = coupling

	for path in b1_paths:
		if types[int(path[0])] == 6 and types[int(path[1])] == 1:
			b1_coupling.append([int(path[0]), int(path[1]), 1, jch_array[int(path[1])]])
		elif types[int(path[1])] == 6 and types[int(path[0])] == 1:
			b1_coupling.append([int(path[0]), int(path[1]), 1, jch_array[int(path[1])]])


	return molid, types, xyz, dist, shift, b1_coupling, b2_coupling, b3_coupling

def jxyz_to_nmrdata(name, path):
	with open(str(path)+str(name), 'r') as f:
		xyzswitch = False
		cplswitch = False
		l = 0

		for line in f:
			l += 1
			items = line.split()

			if l == 1:
				atoms = int(items[0])
				types = np.zeros(atoms, dtype=np.int64)
				xyz = np.zeros((atoms, 3), dtype=np.float64)
				dist = np.zeros((atoms, atoms), dtype=np.int64)
				shift = np.zeros(atoms, dtype=np.float64)
				continue

			if l == 2:
				try:
					molid = items[0]
				except:
					molid = 0

			if 'END' in items:
				xyzswitch = False
				cplswitch = False

			if 'XYZ_DATA' in items:
				xyzswitch = True

			if 'COUPLING_DATA' in items:
				coupling1b = []
				coupling2b = []
				coupling3b = []
				cplswitch = True

			if xyzswitch:
				if len(items) != 5:
					continue
				i = int(items[0])
				types[i] = int(items[1])
				xyz[i][0] = float(items[2])
				xyz[i][1] = float(items[3])
				xyz[i][2] = float(items[4])

			if cplswitch:
				if len(items) < 4:
					continue
				try:
					int(items[1])
				except:
					continue

				if len(items) == 4:
					coupling1b.append([int(items[0]), int(items[1]), int(items[3]), float(items[2])])
				if len(items) == 6:
					coupling2b.append([int(items[1]), int(items[2]), int(items[3]), int(items[4]), float(items[5])])
				if len(items) == 7:
					coupling3b.append([int(items[1]), int(items[2]), int(items[3]), int(items[4]), int(items[5]), float(items[6])])

	b1_coupling = np.asarray(coupling1b)
	b2_coupling = np.asarray(coupling2b)
	b3_coupling = np.asarray(coupling3b)

	return molid, types, xyz, dist, shift, b1_coupling, b2_coupling, b3_coupling

def syg_jxyz_to_nmrdata(name, path):
	with open(str(path)+str(name), 'r') as f:
		xyzswitch = False
		cplswitch = False
		l = 0

		for line in f:
			l += 1
			items = line.split()

			if l == 1:
				atoms = int(items[0])
				types = np.zeros(atoms, dtype=np.int64)
				xyz = np.zeros((atoms, 3), dtype=np.float64)
				dist = np.zeros((atoms, atoms), dtype=np.int64)
				shift = np.zeros(atoms, dtype=np.float64)
				continue

			if l == 2:
				try:
					molid = items[0]
				except:
					molid = 0

			if 'END' in items:
				xyzswitch = False
				cplswitch = False

			if 'XYZ_DATA' in items:
				xyzswitch = True

			if 'COUPLING_DATA' in items:
				coupling1b = []
				coupling2b = []
				coupling3b = []
				cplswitch = True

			if xyzswitch:
				if len(items) != 5:
					continue
				i = int(items[0])
				types[i] = int(items[1])
				xyz[i][0] = float(items[2])
				xyz[i][1] = float(items[3])
				xyz[i][2] = float(items[4])

			if cplswitch:
				if len(items) < 4:
					continue
				try:
					int(items[1])
				except:
					continue

				if len(items) == 6:
					coupling1b.append([int(items[0]), int(items[1]), int(items[4]), float(items[5])])
				if len(items) == 6:
					coupling2b.append([int(items[1]), int(items[2]), int(items[3]), int(items[4]), float(items[5])])
				if len(items) == 7:
					coupling3b.append([int(items[1]), int(items[2]), int(items[3]), int(items[4]), int(items[5]), float(items[6])])

	b1_coupling = np.asarray(coupling1b)
	b2_coupling = np.asarray(coupling2b)
	b3_coupling = np.asarray(coupling3b)

	return molid, types, xyz, dist, shift, b1_coupling, b2_coupling, b3_coupling

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

def raw_to_xyz(file, titleline, xyz, types, append=True):
	periodic_table = Get_periodic_table()
	with open(file, 'w') as f:
		print(len(types), file=f)
		print(titleline, file=f)
		for i in range(len(types)):
			string = "{i:<10s}\t{x:<10.6f}\t{y:<10.6f}\t{z:<10.6f}".format(i=periodic_table[types[i]],
																			x=xyz[i][0],
																			y=xyz[i][1],
																			z=xyz[i][2])
			print(string, file=f)

def pdb_to_nmrdummy(file):

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


	return molid, types, xyz, dist, shift, b1_coupling, b2_coupling, b3_coupling

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


	return molid, types, xyz, dist, shift, b1_coupling, b2_coupling, b3_coupling

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





###

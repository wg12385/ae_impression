import numpy as np
import pybel as pyb
from .pybel_plugins import *
from reference.periodic_table import Get_periodic_table
from reference.tantillo import Get_tantillo_factors
import os

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

def nmrmol_to_pybmol(nmrmol):

	file = nmrmol.molid + 'tmp.xyz'

	nmrmol_to_xyz(nmrmol, file)

	pybmol = next(pyb.readfile('xyz', file))

	os.remove(file)

	return pybmol, nmrmol.molid

def labelmaker(cpl, mol):
	Periodic_table = Get_periodic_table()
	lent = len(cpl)-3
	label = str(lent) + str('J')
	if mol.types[int(cpl[0])] >= mol.types[int(cpl[lent])]:
		label = label + str(Periodic_table[mol.types[int(cpl[0])]]) + str(Periodic_table[mol.types[int(cpl[lent])]])
	else:
		label = label + str(Periodic_table[mol.types[int(cpl[lent])]]) + str(Periodic_table[mol.types[int(cpl[0])]])

	return label

def nmredata_to_nmrmol(file):

	mol = next(pyb.readfile('sdf', file))
	type_list, types = mol_read_type(mol)
	print(types)
	xyz = mol_read_xyz(mol)
	dist = mol_read_dist(mol)

	shift = np.zeros(len(types))
	shift_var = np.zeros(len(types))
	coupling1b = []
	coupling2b = []
	coupling3b = []
	var1b = []
	var2b = []
	var3b = []

	shift_switch = False
	j_switch = False

	ii = 0
	with open(file, 'r') as f:
		for line in f:
			ii += 1
			if ii == 1:
				molid = line.split()[0]


			if 'NMREDATA_ASSIGNMENT' in line:
				shift_switch = True
				j_switch = False

			if 'NMREDATA_J' in line:
				j_switch = True
				shift_switch = False

			if shift_switch:
				items = line.split(',')
				if len(items) != 4:
					continue
				i = int(items[0])
				types[i] = int(items[2])
				shift[i] = float(items[1])
				shift_var[i] = float(items[3].split()[0])

			if j_switch:
				items = line.split(',')
				if len(items) < 5:
					continue
				try:
					int(items[1])
				except:
					continue

				label = items[3].strip()
				length = int(label[0])
				var = items[-1].split()[0]
				if length == 1:
					# (items[0]=label) ||| at1 at2 len coupling
					coupling1b.append([int(items[0].strip()), int(items[1].strip()), length, float(items[2].strip())])
					var1b.append(float(var))

				if length == 2:
					coupling2b.append([int(items[0].strip()), int(items[1].strip()), int(items[2].strip()), length,  float(items[3].strip())])
					var2b.append(float(var))

				if length == 3:
					coupling3b.append([int(items[0].strip()), int(items[1].strip()), int(items[2].strip()), int(items[3].strip()), length,  float(items[4].strip())])
					var3b.append(float(var))

	return molid, types, xyz, dist, shift, shift_var, coupling1b, var1b, coupling2b, var2b, coupling3b, var3b


def nmrmol_to_nmredata(mol, outname):

	pybmol, _ = nmrmol_to_pybmol(mol)
	bond_table = mol_get_bond_table(pybmol)

	periodic_table = Get_periodic_table()
	outfile = outname.split('.')[0] + '.nmredata.sdf'

	lines = []

	bonds = 0
	for at1 in range(len(mol.types)):
		for at2 in range(at1, len(mol.types)):

			if bond_table[at1][at2] >= 1:
				bonds += 1

	atoms = len(mol.types)
	chiral = 0

	lines.append(outname.split('.')[0])
	lines.append('auto-ENRICH')
	lines.append('')

	# Structure section
	string = '{atoms:>3d}{bonds:>3d}  0  0{chiral:>3d}  0  0  0  0  0  1 V2000'.format(atoms=atoms,
																							bonds=bonds,
																							chiral=chiral)
	lines.append(string)

	for i, xyz in enumerate(mol.xyz):
		string = '{x:>10.4f}{y:>10.4f}{z:>10.4f} {typechar:>3s} 0  0  0  0  0  0  0  0  0  0  0  0'.format(x=xyz[0],
																												y=xyz[1],
																												z=xyz[2],
																												typechar = periodic_table[mol.types[i]])
		lines.append(string)

	for at1 in range(len(mol.types)):
		for at2 in range(at1, len(mol.types)):

			if bond_table[at1][at2] >= 1:
				string = '{at1:>3d}{at2:>3d}{bond:>3d}  0  0  0  0'.format(at1=at1+1,
																				at2=at2+1,
																				bond=bond_table[at1][at2])
				lines.append(string)

	lines.append('M\tEND'.format())


	# assignment section
	lines.append('')
	lines.append('> <NMREDATA_ASSIGNMENT>')
	for i, shift, type, var in zip(range(len(mol.types)), mol.shift, mol.types, mol.shift_var):
		string = " {atom:<5d}, {shift:<15.8f}, {type:<5d}, {variance:<15.8f}\\".format(atom=i, shift=shift, type=type, variance=var)
		lines.append(string)

	lines.append('')
	lines.append('> <NMREDATA_J>')

	for i in range(len(mol.coupling1b)):
		if mol.coupling1b[i][0] > mol.coupling1b[i][1]:
			continue
		label = labelmaker(mol.coupling1b[i], mol)
		string = " {a1:<10d}, {a2:<10d}, {coupling:<15.8f}, {label:<10s}, {var:<15.8f}\\".format(a1=int(mol.coupling1b[i][0]),
																	a2=int(mol.coupling1b[i][1]),
																	coupling=mol.coupling1b[i][3],
																	label=label,
																	var=mol.var1b[i])
		lines.append(string)


	for i in range(len(mol.coupling2b)):
		if mol.coupling2b[i][0] > mol.coupling2b[i][2]:
			continue
		label = labelmaker(mol.coupling2b[i], mol)
		string = " {a1:<10d},\t{a2:<10d},\t{a3:<10d},\t{coupling:<15.8f},\t{label:<10s},\t{var:<15.8f}\\".format(label=label,
																				a1=int(mol.coupling2b[i][0]),
																				a2=int(mol.coupling2b[i][1]),
																				a3=int(mol.coupling2b[i][2]),
																				coupling=mol.coupling2b[i][4],
																				var=mol.var2b[i])
		lines.append(string)

	for i in range(len(mol.coupling3b)):
		if mol.coupling3b[i][0] > mol.coupling3b[i][3]:
			continue
		label = labelmaker(mol.coupling3b[i], mol)
		string = " {a1:<10d}, {a2:<10d}, {a3:<10d}, {a4:<10d}, {coupling:<15.8f}, {lbl:<10s}, {var:<15.8f}\\".format(lbl=label,
																				a1=int(mol.coupling3b[i][0]),
																				a2=int(mol.coupling3b[i][1]),
																				a3=int(mol.coupling3b[i][2]),
																				a4=int(mol.coupling3b[i][3]),
																				coupling=mol.coupling3b[i][5],
																				var=mol.var3b[i])
		lines.append(string)

	with open(outfile, 'w') as f:
		for line in lines:
			print(line, file=f)



















####

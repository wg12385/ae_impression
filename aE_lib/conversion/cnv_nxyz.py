# conversions to and from nxyz

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

def nxyz_to_nmrdata(name, path=''):

	with open(str(path)+str(name), 'r') as f:
		xyzswitch = False
		cplswitch = False
		l = 0

		coupling1b = []
		coupling2b = []
		coupling3b = []
		cpl_var1b = []
		cpl_var2b = []
		cpl_var3b = []

		for line in f:
			l += 1
			items = line.split()

			if l == 1:
				atoms = int(items[0])
				types = np.zeros(atoms, dtype=np.int64)
				xyz = np.zeros((atoms, 3), dtype=np.float64)
				dist = np.zeros((atoms, atoms), dtype=np.int64)
				shift = np.zeros(atoms, dtype=np.float64)
				shift_var = np.zeros(atoms, dtype=np.float64)
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
				if len(items) < 6 or len(items) > 7:
					continue
				i = int(items[0])
				types[i] = int(items[1])
				xyz[i][0] = float(items[2])
				xyz[i][1] = float(items[3])
				xyz[i][2] = float(items[4])
				shift[i] = float(items[5])
				if len(items) == 7:
					shift_var[i] = float(items[6])

			if cplswitch:
				if len(items) < 4:
					continue
				try:
					int(items[1])
				except:
					continue

				label = items[0]
				length = int(label[0])

				if length == 1:
					# (items[0]=label) ||| at1 at2 len coupling
					coupling1b.append([int(items[1]), int(items[2]), int(items[3]), float(items[4])])
					try:
						cpl_var1b.append(float(items[5]))
					except:
						cpl_var1b.append(0.0)
				if length == 2:
					coupling2b.append([int(items[1]), int(items[2]), int(items[3]), int(items[4]), float(items[5])])
					try:
						cpl_var2b.append(float(items[6]))
					except:
						cpl_var2b.append(0.0)
				if length == 3:
					coupling3b.append([int(items[1]), int(items[2]), int(items[3]), int(items[4]), int(items[5]), float(items[6])])
					try:
						cpl_var3b.append(float(items[7]))
					except:
						cpl_var3b.append(0.0)

	b1_coupling = np.asarray(coupling1b)
	b2_coupling = np.asarray(coupling2b)
	b3_coupling = np.asarray(coupling3b)

	var1b = np.asarray(cpl_var1b)
	var2b = np.asarray(cpl_var2b)
	var3b = np.asarray(cpl_var3b)

	return molid, types, xyz, dist, shift, shift_var, b1_coupling, var1b, b2_coupling, var2b, b3_coupling, var3b

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
			string = "{i:<10d}\t{typ:<10d}\t{x:<10.6f}\t{y:<10.6f}\t{z:<10.8f}\t{d:<10.8f}\t{v:<10.8f}".format(i=i,
																						typ=int(mol.types[i]),
																						x=mol.xyz[i][0],
																						y=mol.xyz[i][1],
																						z=mol.xyz[i][2],
																						d=mol.shift[i],
																						v=mol.shift_var[i])
			print(string, file=f)
		print('END', file=f)
		print('', file=f)
		print('COUPLING_DATA', file=f)

		for i in range(len(mol.coupling1b)):
			if mol.coupling1b[i][0] > mol.coupling1b[i][1]:
				continue
			label = labelmaker(mol.coupling1b[i], mol)
			string = "{lbl:<10s}\t{a1:<10d}\t{a2:<10d}\t{len:<10d}\t{coupling:<10.8f}\t{var:<10.8f}".format(lbl=label,
																					a1=int(mol.coupling1b[i][0]),
																					a2=int(mol.coupling1b[i][1]),
																					len=int(mol.coupling1b[i][2]),
																					coupling=mol.coupling1b[i][3],
																					var=mol.var1b[i])
			print(string, file=f)

		for i in range(len(mol.coupling2b)):
			if mol.coupling2b[i][0] > mol.coupling2b[i][2]:
				continue
			label = labelmaker(mol.coupling2b[i], mol)
			string = "{lbl:<10s}\t{a1:<10d}\t{a2:<10d}\t{a3:<10d}\t{len:<10d}\t{coupling:<10.8f}\t{var:<10.8f}".format(lbl=label,
																					a1=int(mol.coupling2b[i][0]),
																					a2=int(mol.coupling2b[i][1]),
																					a3=int(mol.coupling2b[i][2]),
																					len=int(mol.coupling2b[i][3]),
																					coupling=mol.coupling2b[i][4],
																					var=mol.var2b[i])
			print(string, file=f)

		for i in range(len(mol.coupling3b)):
			if mol.coupling3b[i][0] > mol.coupling3b[i][3]:
				continue
			label = labelmaker(mol.coupling3b[i], mol)
			string = "{lbl:<10s}\t{a1:<10d}\t{a2:<10d}\t{a3:<10d}\t{a4:<10d}\t{len:<10d}\t{coupling:<10.8f}\t{var:<10.8f}".format(lbl=label,
																					a1=int(mol.coupling3b[i][0]),
																					a2=int(mol.coupling3b[i][1]),
																					a3=int(mol.coupling3b[i][2]),
																					a4=int(mol.coupling3b[i][3]),
																					len=int(mol.coupling3b[i][4]),
																					coupling=mol.coupling3b[i][5],
																					var=mol.var3b[i])
			print(string, file=f)
		print('END', file=f)

	return outname

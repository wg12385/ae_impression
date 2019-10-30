
from util.labelmaker import labelmaker
from reference.periodic_table import Get_periodic_table

def nmrmol_to_nmredata(mol, outname):

	pybmol, _ = nmrmol_to_pybmol(mol)
	bond_table = mol_get_bond_table(pybmol)

	periodic_table = Get_periodic_table()

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

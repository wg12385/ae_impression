from reference.periodic_table import Get_periodic_table

def labelmaker(i, j, mol):
	Periodic_table = Get_periodic_table()
	lent = mol.coupling_len[i][j]
	label = str(lent) + str('J')
	if mol.types[int(i)] >= mol.types[int(j)]:
		label = label + str(Periodic_table[mol.types[int(i)]]) + str(Periodic_table[mol.types[int(j)]])
	else:
		label = label + str(Periodic_table[mol.types[int(j)]]) + str(Periodic_table[mol.types[int(i)]])

	return label

from reference.periodic_table import Get_periodic_table

def labelmaker(cpl, mol):
	Periodic_table = Get_periodic_table()
	lent = len(cpl)-3
	label = str(lent) + str('J')
	if mol.types[int(cpl[0])] >= mol.types[int(cpl[lent])]:
		label = label + str(Periodic_table[mol.types[int(cpl[0])]]) + str(Periodic_table[mol.types[int(cpl[lent])]])
	else:
		label = label + str(Periodic_table[mol.types[int(cpl[lent])]]) + str(Periodic_table[mol.types[int(cpl[0])]])

	return label

from reference.periodic_table import Get_periodic_table

def flag_to_target(flag):

	p_table = Get_periodic_table()
	if len(label) == 4:
		if str(label[0]) == 'J':
			length = int(label[1])
		else:
			length = int(label[0])
		atype1 = int(p_table.index(str(label[2])))
		atype2 = int(p_table.index(str(label[3])))

		if atype1 >= atype2:
			return [length, atype1, atype2]
		else:
			return [length, atype2, atype1]

	elif len(label) == 3:
		atype = int(p_table.index(str(label[0])))
		return [atype]

	else:
		print('label, ', label, ' not recognised, coupling label format is <nJxy> . . .')
		print('label, ', label, ' not recognised, chemical shift label format is <XCS> . . .')
		return 0

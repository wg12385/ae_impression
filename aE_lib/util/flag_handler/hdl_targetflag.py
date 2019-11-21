from reference.periodic_table import Get_periodic_table

def flag_to_target(flag):

	p_table = Get_periodic_table()
	if len(flag) == 4:
		if str(flag[0]) == 'J':
			length = int(flag[1])
		else:
			length = int(flag[0])
		atype1 = int(p_table.index(str(flag[2])))
		atype2 = int(p_table.index(str(flag[3])))

		if atype1 >= atype2:
			return [length, atype1, atype2]
		else:
			return [length, atype2, atype1]

	elif len(flag) == 3:
		atype = int(p_table.index(str(flag[0])))
		return [atype]

	else:
		print('flag, ', flag, ' not recognised, coupling flag format is <nJxy> . . .')
		print('flag, ', flag, ' not recognised, chemical shift flag format is <XCS> . . .')
		return 0
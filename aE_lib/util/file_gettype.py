


def get_type(filename):

	extension = filename.split('.')[-1]
	if extension == 'sdf':
		if len(filename.split('.')) > 2:
			if filename.split('.')[-2] == 'nmredata':
				type = 'nmredata'
			else:
				type = 'sdf'
		else:
			type = 'sdf'

	elif extension == 'xyz':
		type = 'xyz'

	elif extension == 'log':
		type = 'g09'

	elif extension == 'mol2':
		type = 'mol2'

	else:
		print('type not recognised for file ', filename, ' please check file and specify type')
		type = ''

	return type

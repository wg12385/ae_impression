
def print_file(filename, strings):
	with open(filename, 'w') as f:
		for string in strings:
			print(string, file=f)



def get_unique_part(files):

	parts = []

	for file in files:
		parts.append(file.split('/')[-1].split('.')[0].split('_'))

	label_part = 0

	for p1, part1 in enumerate(parts):
		for p2, part2 in enumerate(parts):
			if p1 == p2:
				continue
			for i in range(len(part1)):
				if part1[i] != part2[i]:
					label_part = i


	return label_part

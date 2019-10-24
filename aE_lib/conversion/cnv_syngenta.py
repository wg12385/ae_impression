# conversion for data provided by syngenta



def sygcml_to_nmrdata(name, path):

	mol = next(pyb.readfile('cml', str(path+name)))
	atoms = len(mol.atoms)
	xyz = mol_read_xyz(mol)
	type_list, types = mol_read_type(mol)
	shift = np.zeros(atoms, dtype=np.float64)
	dist = np.zeros((atoms,atoms), dtype=np.float64)

	b1_paths = []
	b2_paths = []
	b3_paths = []

	b1_coupling = []
	b2_coupling = []
	b3_coupling = []

	for x in set(types):
		for y in set(types):
			b1_paths.extend(mol_pathway_finder(mol, x, y, 1))
			b2_paths.extend(mol_pathway_finder(mol, x, y, 2))
			b3_paths.extend(mol_pathway_finder(mol, x, y, 3))

	j_array = np.zeros((atoms, atoms), dtype=np.float64)

	jch_array = np.zeros(atoms, dtype=np.float64)

	with open(str(path+name), 'r') as f:
		for line in f:
			items = line.split()
			if '<molecule id=' in line:
				molid = int(line.split('<molecule id="')[1].split('_')[0])
			if '<atom id=' in line:
				atomname = line.split('id="')[1].split('"')[0]
				type = line.split('elementType="')[1].split('"')[0]
				x = float(line.split('x3="')[1].split('"')[0])
				y = float(line.split('y3="')[1].split('"')[0])
				z = float(line.split('z3="')[1].split('"')[0])
			if '<scalar title="Exp1JCH"' in line:
			#if '<scalar title="NWCHEM_1JCH"' in line:
				try:
					coupling = float(line.split('dataType="xsd:string">')[1].split('</scalar>')[0])
				except Exception as e:
					print(name, e, line)
					continue

				for i in range(atoms):
					if x == xyz[i][0] and y == xyz[i][1] and z == xyz[i][2]:
						jch_array[i] = coupling

	for path in b1_paths:
		if types[int(path[0])] == 6 and types[int(path[1])] == 1:
			b1_coupling.append([int(path[0]), int(path[1]), 1, jch_array[int(path[1])]])
		elif types[int(path[1])] == 6 and types[int(path[0])] == 1:
			b1_coupling.append([int(path[0]), int(path[1]), 1, jch_array[int(path[1])]])

	shift_var = np.zeros(len(shift), dtype=np.float64)
	var1b = np.zeros(len(b1_coupling), dtype=np.float64)
	var2b = np.zeros(len(b2_coupling), dtype=np.float64)
	var3b = np.zeros(len(b3_coupling), dtype=np.float64)

	return molid, types, xyz, dist, shift, shift_var, b1_coupling, var1b, b2_coupling, var2b, b3_coupling, var3b

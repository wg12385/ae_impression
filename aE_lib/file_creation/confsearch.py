
# Write file to perform conformational search
def make_confsearch_script(scriptname, pickle_file, molecule, smiles, aE_dir='../../aE_lib/', path='', iterations=2000,
							RMSthresh=10.0, maxconfs=10, Ethresh=999.99):

	strings = []

	strings.append('import sys')
	strings.append("sys.path.append('{dir:<10s}')".format(dir=aE_dir))
	strings.append('import pickle')
	strings.append('from file_creation.structure_formats import xyz')
	strings.append("file = '{0:<10s}'".format(pickle_file.split('/')[-1]))
	strings.append("molecule = pickle.load(open(file, 'rb'))")
	strings.append("molecule.generate_conformers('{smiles:<10s}', path='{path:<10s}', iterations={its:<10d}".format(smiles=smiles,
																												path=path,
																												its=iterations))
	strings.append(", RMSthresh={RMS:<10f}, maxconfs={maxconfs:<10d}, Ethresh={Ethresh:<10f})".format(RMS=RMSthresh,
																									maxconfs=maxconfs,
																									Ethresh=Ethresh))


	strings.append("for conformer in molecule.conformers:")
	strings.append("\txyz_file = 'conf_search/' +  conformer.molid + '.xyz'")
	strings.append("\txyz.nmrmol_to_xyz(molecule, xyz_file)")
	strings.append("pickle.dump(molecule, open(file, 'wb'))")


	with open(path + scriptname, 'w') as f:
		for string in strings:
			print(string, file=f)

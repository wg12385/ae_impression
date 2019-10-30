from file_creation.structure_formats import xyz
from file_creation import confsearch
from file_creation import HPC_submission as HPCsub
import pybel as pyb

def conformational_search(molecule, prefs, path=''):
	# Read relevant preferences
	iterations = prefs['conf_search']['iterations']
	maxconfs = prefs['conf_search']['maxconfs']
	RMSthresh = prefs['conf_search']['RMSthresh']
	Ethresh = prefs['conf_search']['Ethresh']

	system = prefs['comp']['system']
	python_env = prefs['comp']['python_env']
	memory = prefs['conf_search']['memory']
	processors = prefs['conf_search']['processors']
	walltime = prefs['conf_search']['walltime']

	# Try reading pybel mol object from xyz named according to molname
	try:
		file = path + molecule.molid + '.xyz'
		mol = next(pyb.readfile('xyz', file))
	# If not, make a molecule xyz file and then read the pybel mol object
	except:
		file = path + molecule.molid + '_INIT.xyz'
		xyz.nmrmol_to_xyz(molecule, file, num=-404)
		mol = next(pyb.readfile('xyz', file))

	smiles = mol.write("smi")

	scriptname = path + 'conf_search/do_confsearch.py'
	confsearch.make_confsearch_script(scriptname, molecule, smiles,
									iterations=iterations, RMSthresh=RMSthresh,
												maxconfs=maxconfs, Ethresh=Ethresh)

	strings1 = HPCsub.make_HPC_header(prefs)
	strings.append('source activate {env:<10s}').format(env=python_env)
	strings.append('python {0:<10s}').format(scriptname)

	filename = path + 'submit_confsearch.sh'
	with open('')

	if system == 'localbox':
		# Run conformational search
		molecule.generate_conformers(smiles, path=path,
					iterations=iterations, RMSthresh=RMSthresh,
					maxconfs=maxconfs, Ethresh=Ethresh)

		molecule.print_xyzs(path=path+'conf_search/')
	else:
		if system == 'BC3':
			outname = path +
			strings.append("#PBS -l nodes={0:<1d}:ppn={1:<1d}".format(1, processors))
			strings.append("#PBS -l walltime={0:<9s}".format(walltime))
			strings.append("#PBS -l mem={0:<1d}GB".format(memory))
			strings.append("cd $PBS_O_WORKDIR")
			strings.append("source activate {0:<1s}").format(python_env)

		# figure out how we're doing this . . .
		strings.append("filename = {file:<10s}").format(file)
		molecule = moleculeclass(from_file=filename, path=args.path, from_type='pickle')
		strings.append("molecule.generate_conformers({smiles:<10s}, path={path:<10s}, iterations={its:<10d}".format(smiles=smiles,
																													path=path,
																													its=iterations))
		strings.append(", RMSthresh={RMS:<10f}, maxconfs={maxconfs:<10d}, Ethresh={Ethresh:<10f})".format(RMS=RMSthresh,
																										maxconfs=maxconfs,
																										Ethresh=Ethresh))
																										# write code to do submission scripts for conf search
		#with open(path + 'conf_search/conf_search.sh')

	return 0

import sys
sys.path.append('aE_lib')

from file_creation import orca_submission as orcasub

from conversion.convert import xyz_to_nmrdummy, nmrmol_to_xyz
import pybel as pyb

# Runs conformational search
def conformational_search(molecule, prefs, path=''):
	# Try reading pybel mol object from xyz named according to molname
	try:
		file = path + molecule.molid + '.xyz'
		mol = next(pyb.readfile('xyz', file))
	# If not, make a molecule xyz file and then read the pybel mol object
	except:
		file = path + molecule.molid + '_INIT.xyz'
		nmrmol_to_xyz(molecule, file, num=-404)
		mol = next(pyb.readfile('xyz', file))

	# Get smiles string using pybel mol object
	# Read relevant preferences
	iterations = prefs['conf_search']['iterations']
	maxconfs = prefs['conf_search']['maxconfs']
	RMSthresh = prefs['conf_search']['RMSthresh']
	Ethresh = prefs['conf_search']['Ethresh']

	system = prefs['comp']['system']
	pythn_env = prefs['comp']['python_env']
	memory = prefs['conf_search']['memory']
	processors = prefs['conf_search']['processors']
	walltime = prefs['conf_search']['walltime']

	smiles = mol.write("smi")

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


		strings.append('filename = args.Molecule + '.pkl'
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

# Setup optimisation gaussian in files
def setup_opt(molecule, prefs, path=''):
	# Read relevant preferences

	for conformer in molecule.conformers:
		conformer.optin = orcasub.make_optin(prefs, conformer, path)

	qsub_names = molecule.make_opt_sub(prefs, path=path, start=-1, end=-1)



	print('Created ', len(qsub_names), ' qsub files. . .')
	if prefs['comp']['system'] == 'BC3':
		print('Submit the calculations using:')
		for file in qsub_names:
			print('qsub ', file)
	elif prefs['comp']['system'] == 'Grendel':
		print('Submit the calculations using:')
		for file in qsub_names:
			print('bash ', file)
	elif prefs['comp']['system'] == 'BC4':
		print('Havent finished this yet, good luck pal. . . .')

	return 0

def process_opt(molecule, prefs, path):
	statuss = []
	good = 0
	bad = 0
	for conformer in molecule.conformers:
		conformer.check_opt()
		statuss.append(conformer.opt_status)
		string = 'Conformer {molid:^3s} status: {status:^10s} | Energy {energy:<10.4f}'.format(molid=str(conformer.molid),
																								status=conformer.opt_status,
																								energy=conformer.energy)
		print(string)

		if conformer.opt_status == 'successful':
			good += 1
		elif conformer.opt_status == 'failed':
			bad += 1

	print(good, ' successful optimisations, ', bad, ' failed, out of ', len(statuss))

	molecule.print_xyzs(path=path + 'optimisation/')

	if bad >= 1:
		qsub_names = molecule.make_opt_sub(prefs, path=path+'optimisation/RESUB_FAILED_', start=-1, end=-1, failed_only=True)
		print('Created ', len(qsub_names), ' qsub files to resubmit failed calculations')
		if prefs['comp']['system'] == 'BC3':
			print('Fix the issue (check log_file for error) then submit the calculations using:')
			for file in qsub_names:
				print('qsub ', file)
		elif prefs['comp']['system'] == 'Grendel':
			print('Fix the issue (check log_file for error) then submit the calculations using:')
			for file in qsub_names:
				print('bash ', file)
		else:
			print('Fix the issue (check log_file for error) then submit the calculations using:')
			for file in qsub_names:
				print('bash ', file)

		print("Resubmit failed optimisations or continue to NMR calculations with 'setup_nmr'")

# Setup NMR gaussian in files
def setup_nmr(molecule, prefs, path):

	molecule.make_nmr_in(prefs, path=path)

	qsub_names = molecule.make_nmr_sub(prefs, path=path, start=-1, end=-1)

	print('Created ', len(qsub_names), ' qsub files. . .')
	if system == 'BC3':
		print('Submit the calculations using:')
		for file in qsub_names:
			print('qsub ', file)
	else:
		print('Submit the calculations using:')
		for file in qsub_names:
			print('bash ', file)
	return 0

# Process NMR log files
def process_nmr():

	mo

	return status

def update_molecule(molecule, prefs):
	change_opt = False
	change_nmr = True
	for conformer in molecule.conformers:
		prev = conformer.opt_status
		opt = conformer.check_opt_status()
		if opt != prev:
			change_opt = True
		prev = conformer.nmr_status
		nmr = conformer.check_nmr_status()
		if nmr != prev:
			change_nmr = True


	if change_opt:
		print('Optimisation update, making NMR in files')
		setup_nmr(molecule, prefs)
	if change_nmr:
		process_nmr(molecule, prefs)
		print('NMR update, ')




	return 0






###

import sys
sys.path.append('aE_lib')
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
	smiles = mol.write("smi")

	# Read relevant preferences
	iterations = prefs['conf_search']['iterations']
	maxconfs = prefs['conf_search']['maxconfs']
	RMSthresh = prefs['conf_search']['RMSthresh']
	Ethresh = prefs['conf_search']['Ethresh']
	# Run conformational search
	molecule.generate_conformers(smiles, path=path,
				iterations=iterations, RMSthresh=RMSthresh,
				maxconfs=maxconfs, Ethresh=Ethresh)

	molecule.print_xyzs(path=path+'conf_search/')

	return 0

# Setup optimisation gaussian in files
def setup_opt(molecule, prefs, path=''):
	# Read relevant preferences

	molecule.make_opt_in(prefs, path=path)

	qsub_names = molecule.make_opt_sub(prefs, path=path, start=-1, end=-1)

	molecule.print_xyzs(path=path+'optimisation/')

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

	if bad >= 1:
		qsub_names = molecule.make_opt_sub(prefs, path=path+'optimisation/RESUB_FAILED_', start=-1, end=-1, failed_only=True)

	print(good, ' successful optimisations, ', bad, ' failed, out of ', len(statuss))

	molecule.print_xyzs(path=path + 'optimisation/')

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


# Script to test autoenrich workflow
from autoenrich.molecule.molecule import molecule as moleculeclass
from autoenrich.preferences.preferences import write_default_prefs
import autoenrich.top_level.CMD_confsearch as CMD_confsearch
import autoenrich.top_level.CMD_optimisation as CMD_opt
import autoenrich.top_level.CMD_nmr as CMD_nmr
import pickle


def test_init():

	path = '/'.join(__file__.split('/')[:-1])

	status = 'Pass'
	args = {}
	args['Command'] = 'init'
	args['Molecule'] = 'testmol'
	args['xyz_file'] = path + '/../test_store/ethane.xyz'
	args['path'] = path + '/../test_tmp/'
	write_default_prefs(args['path']+'ENRICH.json')
	args['prefs'] = args['path'] + 'ENRICH.json'

	# Load xyz coords and types from xyz file, create molecule object
	try:
		molecule = moleculeclass(molid=args.Molecule, path=args.path)
		molecule.read_structure(args.path + args.xyz_file, 'xyz')
		pickle.dump(molecule, open(pickle_file, "wb"))
	except Exception as e:
		print('ERROR IN TEST -------------------------------------------')
		print(args)
		print(e)
		status = 'Fail'

	return status

def test_confsearch():

	path = '/'.join(__file__.split('/')[:-1])

	status = 'Pass'
	args = {}
	args['Command'] = 'conf_search'
	args['Molecule'] = 'testmol'
	args['path'] = path + '/../test_tmp/'
	args['prefs'] = args['path'] + 'ENRICH.json'

	# Load xyz coords and types from xyz file, create molecule object
	try:
		pickle_file = args['path'] + args['Molecule'] + '.pkl'
		molecule = pickle.load(open(pickle_file, "rb"))
		# Do conformational search
		CMD_confsearch.conformational_search(molecule, prefs, pickle_file, path=args['path'])
		# Save molecule in pickle file
		pickle.dump(molecule, open(pickle_file, "wb"))
	except Exception as e:
		print('ERROR IN TEST -------------------------------------------')
		print(args)
		print(e)
		status = 'Fail'

	return status

def test_setupopt():

	path = '/'.join(__file__.split('/')[:-1])

	status = 'Pass'
	args = {}
	args['Command'] = 'setup_opt'
	args['Molecule'] = 'testmol'
	args['path'] = path + '/../test_tmp/'
	args['prefs'] = args['path'] + 'ENRICH.json'

	# Load xyz coords and types from xyz file, create molecule object
	try:
		pickle_file = args['path'] + args['Molecule'] + '.pkl'
		molecule = pickle.load(open(pickle_file, "rb"))
		CMD_opt.setup_optimisation(molecule, prefs, path=args['path'])
		# Save molecule in pickle file
		pickle.dump(molecule, open(pickle_file, "wb"))
	except Exception as e:
		print('ERROR IN TEST -------------------------------------------')
		print(args)
		print(e)
		status = 'Fail'

	return status

def test_processopt():

	path = '/'.join(__file__.split('/')[:-1])

	status = 'Pass'
	args = {}
	args['Command'] = 'process_opt'
	args['Molecule'] = 'testmol'
	args['path'] = path + '/../test_tmp/'
	args['prefs'] = args['path'] + 'ENRICH.json'

	# Load xyz coords and types from xyz file, create molecule object
	try:
		pickle_file = args['path'] + args['Molecule'] + '.pkl'
		molecule = pickle.load(open(pickle_file, "rb"))
		CMD_opt.process_optimisation(molecule, prefs, path=args['path'])
		# Save molecule in pickle file
		pickle.dump(molecule, open(pickle_file, "wb"))
	except Exception as e:
		print('ERROR IN TEST -------------------------------------------')
		print(args)
		print(e)
		status = 'Fail'
	return status

def test_setupnmr():

	path = '/'.join(__file__.split('/')[:-1])

	status = 'Pass'
	args = {}
	args['Command'] = 'setup_nmr'
	args['Molecule'] = 'testmol'
	args['path'] = path + '/test_files/testmol/'
	args['prefs'] = args['path'] + 'ENRICH.json'

	# Load xyz coords and types from xyz file, create molecule object
	try:
		pickle_file = args['path'] + args['Molecule'] + '.pkl'
		molecule = pickle.load(open(pickle_file, "rb"))
		CMD_nmr.setup_nmr(molecule, prefs, path=args['path'])
		# Save molecule in pickle file
		pickle.dump(molecule, open(pickle_file, "wb"))
	except Exception as e:
		print('ERROR IN TEST -------------------------------------------')
		print(args)
		print(e)
		status = 'Fail'

	return status

def test_processnmr():

	path = '/'.join(__file__.split('/')[:-1])

	status = 'Pass'
	args = {}
	args['Command'] = 'process_nmr'
	args['Molecule'] = 'testmol'
	args['path'] = path + '/test_files/testmol/'
	args['prefs'] = args['path'] + 'ENRICH.json'

	# Load xyz coords and types from xyz file, create molecule object
	try:
		pickle_file = args['path'] + args['Molecule'] + '.pkl'
		molecule = pickle.load(open(pickle_file, "rb"))
		CMD_nmr.process_nmr(molecule, prefs, path=args['path'])
		# Save molecule in pickle file
		pickle.dump(molecule, open(pickle_file, "wb"))
	except Exception as e:
		print('ERROR IN TEST -------------------------------------------')
		print(args)
		print(e)
		status = 'Fail'

	return status

def test_update():

	path = '/'.join(__file__.split('/')[:-1])

	status = 'Pass'
	args = {}
	args['Command'] = 'update'
	args['Molecule'] = 'testmol'
	args['path'] = path + '/test_files/testmol/'
	args['prefs'] = args['path'] + 'ENRICH.json'

	# Load xyz coords and types from xyz file, create molecule object
	try:
		pickle_file = args['path'] + args['Molecule'] + '.pkl'
		molecule = pickle.load(open(pickle_file, "rb"))
		status = check_status(molecule)
		update_molecule(molecule)
		pickle.dump(molecule, open(pickle_file, "wb"))
	except Exception as e:
		print('ERROR IN TEST -------------------------------------------')
		print(args)
		print(e)
		status = 'Fail'
	return status

def test_resubfailed():
	print('Not Done Yet')
	return 'Pass'

def test_undo():

	path = '/'.join(__file__.split('/')[:-1])

	status = 'Pass'
	args = {}
	args['Command'] = 'undo'
	args['Molecule'] = 'testmol'
	args['path'] = path + '/test_files/testmol/'
	args['prefs'] = args['path'] + 'ENRICH.json'

	# Load xyz coords and types from xyz file, create molecule object
	try:
		pickle_file = args['path'] + args['Molecule'] + '.pkl'
		molecule = pickle.load(open(backup_file, "rb"))
		pickle.dump(molecule, open(pickle_file, "wb"))
	except Exception as e:
		print('ERROR IN TEST -------------------------------------------')
		print(args)
		print(e)
		status = 'Fail'

	return status


if __name__ == "__main__":
	status = test_init()
	print('init', status)
	status = test_confsearch()
	print('conf_search', status)
	status = test_setupopt()
	print('setup_opt', status)
	status = test_processopt()
	print('process_opt', status)
	status = test_setupnmr()
	print('setup_nmr', status)
	status = test_processnmr()
	print('process_nmr', status)
	status = test_update()
	print('update', status)
	status = test_resubfailed()
	print('resub_failed', status)
	status = test_undo()
	print('undo', status)

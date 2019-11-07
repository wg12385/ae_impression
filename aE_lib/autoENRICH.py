# Central command file for running auto-ENRICH
# command line options run the code
import argparse
#import sys
#sys.path.append('aE_lib')
from preferences.preferences import read_prefs, write_default_prefs
from molecule.molecule import molecule as moleculeclass
from util.header import print_header_aE
import pybel as pyb
import os
import glob

import pickle

import top_level.CMD_confsearch as CMD_confsearch
import top_level.CMD_optimisation as CMD_opt
import top_level.CMD_nmr as CMD_nmr
import top_level.CMD_resubmission as CMD_resub
from top_level.progress import check_proceed
from top_level.user_util import yes_or_no

# Main auto-ENRICH program
if __name__ == "__main__":
	# Argparser
	parser = argparse.ArgumentParser(description='auto-ENRICH')
	# Define name of molecule (used for save/load of pickle file)
	parser.add_argument('Molecule', help='name of molecule')
	# Define command to perform
	parser.add_argument('Command', help='auto-ENRICH command', choices=['undo', 'init', 'conf_search', 'setup_opt', 'process_opt', 'setup_nmr', 'process_nmr', 'update', 'check_status'])
	# Optional arguments
	parser.add_argument('--xyz', help='xyz file to initialise molecule', action="store", dest='xyz_file', default='None')
	parser.add_argument('--path', help='path to molecule pickle file', action="store", dest='path', default='')
	parser.add_argument('--prefs', help='preferences file to use', action="store", dest='prefs', default='ENRICH.json')
	# Parse arguments into args object
	args = parser.parse_args()
	# Print raw input arguments to user
	print(args.Molecule, args.Command)
	# Print pretty banner
	print_header_aE()

	# make sure path is directory
	if args.path[-1] != "/":
		args.path = args.path + "/"

	pickle_file = args.path + args.Molecule + '.pkl'
	backup_file = args.path + args.Molecule + 'BACKUP.pkl'

	# Check for / Read preferences file
	pref_file = args.path+args.prefs
	if os.path.isfile(pref_file):
		prefs = read_prefs(pref_file)
	else:
		# If no preferences file found, warn user, print default file then quit
		print('Could not find preferences file located at: ', pref_file)
		print('Creating default preferences file in ', args.path+'ENRICH.json')
		write_default_prefs(args.path+'ENRICH.json')
		print('Edit preferences file then run again')
		sys.exit(0)

	# Initialisation command
	if args.Command == 'init':
		# Check for supplied xyz file, quit if not found
		if args.xyz_file == 'None':
			print('ERROR: Must supply xyz file to initialise molecule')
			sys.exit(0)

		if os.path.isfile(pickle_file):
			print('Molecule', args.Molecule, 'already exists, do you want to overwrite this molecule ?')
			answer = yes_or_no()
			if not answer:
				sys.exit(0)
			else:
				print('Overwriting molecule. . .')

		# Load xyz coords and types from xyz file, create molecule object
		molecule = moleculeclass(molid=args.Molecule, path=args.path)
		molecule.read_structure(args.path + args.xyz_file, 'xyz')
		pickle.dump(molecule, open(pickle_file, "wb"))


		#molecule = moleculeclass(init_xyz, init_types, name=args.Molecule, path=args.path)
	# If not initialising, get molecule from file
	else:
		# Load molecule object
		molecule = pickle.load(open(pickle_file, "rb"))


	# Get molecule status, print to user
	status = molecule.stage
	print('Molecule ', molecule.molid, ' stage: ', status)
	'''
	proceed, molecule = check_proceed(status, args.Command, molecule, args.path)
	if not proceed:
		print('Exiting. . .')
		sys.exit(0)
	'''

	if args.Command == 'undo':
		molecule = pickle.load(open(backup_file, "rb"))
		pickle.dump(molecule, open(pickle_file, "wb"))
	else:
		pickle.dump(molecule, open( backup_file, "wb"))

	# Conformational Search command
	if args.Command == 'conf_search':
		# If molecule is at the wrong stage, give option to overwrite data

		print('Running conformational search for molecule ', args.Molecule)
		# Make conf_search directory
		try:
			os.mkdir(args.path+'conf_search')
		except:
			pass

		# Do conformational search
		CMD_confsearch.conformational_search(molecule, prefs, pickle_file, path=args.path)
		# Save molecule in pickle file
		pickle.dump(molecule, open(pickle_file, "wb"))


	elif args.Command == 'setup_opt':
		print('Generating optimisation ORCA input files for molecule, ', args.Molecule)
		try:
			os.mkdir(args.path+'optimisation/')
		except:
			pass
		CMD_opt.setup_optimisation(molecule, prefs, path=args.path)
		molecule.stage = 'opt'
		pickle.dump(molecule, open(pickle_file, "wb"))

	elif args.Command == 'process_opt':
		print('Processing optimisation log files')
		CMD_opt.process_optimisation(molecule, prefs, path=args.path)

		pickle.dump(molecule, open(pickle_file, "wb"))

	elif args.Command == 'setup_nmr':
		print('Generating NMR ORCA input files for molecule, ', args.Molecule)
		try:
			os.mkdir(args.path+'NMR')
		except:
			pass
		CMD_nmr.setup_nmr(molecule, prefs, path=args.path)
		molecule.stage = 'nmr'
		pickle.dump(molecule, open(pickle_file, "wb"))

	elif args.Command == 'process_nmr':
		try:
			os.mkdir(args.path+'OUTPUT/')
		except:
			pass
		CMD_nmr.process_nmr(molecule)
		molecule.stage = 'post'
		pickle.dump(molecule, open(pickle_file, "wb"))


	elif args.Command == 'resub_failed':
		print('Preparing files for resubmission')
		CMD_resub.setup_resubmission(molecule, prefs, path=args.path)
		pickle.dump(molecule, open(pickle_file, "wb"))


	elif args.Command == 'check_status':
		pickle.dump(molecule, open(pickle_file, "wb"))

	elif args.Command == 'update':
		status = check_status(molecule)
		update_molecule(molecule)
		pickle.dump(molecule, open(pickle_file, "wb"))
















###

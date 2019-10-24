# Central command file for running auto-ENRICH
# command line options run the code
import argparse
import sys
sys.path.append('aE_lib')
from preferences.preferences import read_prefs, write_default_prefs
from molecule.molecule import molecule as moleculeclass
from util.header import print_header
import pybel as pyb
import os
import glob

from top_level.commands import *
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
	print_header()
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

		filename = args.path + args.Molecule + '.pkl'
		if os.path.isfile(filename):
			print('Molecule', args.Molecule, 'already exists, do you want to overwrite this molecule ?')
			answer = yes_or_no()
			if not answer:
				sys.exit(0)
			else:
				print('Overwriting molecule. . .')

		# Load xyz coords and types from xyz file, create molecule object
		molecule = moleculeclass(path=args.path, from_file=args.xyz_file, from_type='xyz')
		molecule.save_molecule_to_file()


		#molecule = moleculeclass(init_xyz, init_types, name=args.Molecule, path=args.path)
	# If not initialising, get molecule from file
	else:
		# Load molecule object
		filename = args.Molecule + '.pkl'
		molecule = moleculeclass(from_file=filename, path=args.path, from_type='pickle')

	# Get molecule status, print to user
	status = molecule.stage
	print('Molecule ', molecule.molid, ' stage: ', status)
	proceed, molecule = check_proceed(status, args.Command, molecule, args.path)
	if not proceed:
		print('Exiting. . .')
		sys.exit(0)

	if args.Command == 'undo':
		molecule.load_molecule(path=args.path+'BACKUP_')
	else:
		molecule.save_molecule_to_file(path=args.path+'BACKUP_')

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
		conformational_search(molecule, prefs, path=args.path)
		# Save molecule in pickle file
		molecule.save_molecule_to_file(path=args.path)
		print('Conformational search complete, number of conformers generated:', len(molecule.conformers))


	elif args.Command == 'setup_opt':
		print('Generating optimisation ORCA input files for molecule, ', args.Molecule)
		try:
			os.mkdir(args.path+'optimisation')
		except:
			pass
		setup_opt(molecule, prefs, path=args.path)
		molecule.save_molecule_to_file(path=args.path)

	elif args.Command == 'process_opt':
		print('Processing optimisation log files')
		process_opt(molecule, prefs, path=args.path)

		molecule.save_molecule_to_file(path=args.path)


	elif args.Command == 'setup_nmr':
		print('Generating NMR ORCA input files for molecule, ', args.Molecule)
		try:
			os.mkdir(args.path+'NMR')
		except:
			pass
		setup_nmr(molecule, prefs, path=args.path)
		molecule.save_molecule_to_file(path=args.path)

	elif args.Command == 'process_nmr':
		process_nmr(molecule)
		molecule.save_molecule_to_file(path=args.path)

	elif args.Command == 'check_status':
		molecule.save_molecule_to_file(path=args.path)

	elif args.Command == 'update':
		status = check_status(molecule)
		update_molecule(molecule)
		molecule.save_molecule_to_file(path=args.path)

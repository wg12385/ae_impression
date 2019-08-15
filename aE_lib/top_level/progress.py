from top_level.user_util import yes_or_no

def check_proceed(stage, command, molecule):
	proceed = False

	if command == 'init':
		proceed = True

	elif command == 'conf_search':
		if status != 'init':
			print('Conformational search already run for molecule,', args.Molecule, '. . . continuing will erase existing molecule data.')
			print('Do you want to continue ?')
			# Get user conformation for overwrite
			answer = yes_or_no()
			if answer:
				# enumerate files in conf_search directory then delete them
				files_to_remove = glob.glob(args.path + 'conf_search/*')
				for file in files_to_remove:
					os.remove(file)
				# remove conf_search directory
				os.rmdir(args.path + 'conf_search')
				# Re-initialise molecule from existing molecule, retaining only xyz and types
				init_xyz = molecule.xyz
				init_types = molecule.types
				molecule = moleculeclass(init_xyz, init_types, name=args.Molecule, path=args.path)
				proceed = True
			else:
				proceed = False

	elif command == 'setup_opt':
		if stage in ['init', 'pre-opt']:
			proceed = True
		else:
			print('Already setup optimisations, do you want to continue ? This may overwrite/remove existing com/log files')
			answer = yes_or_no()
			if answer:
				molecule.stage='pre-opt'
				proceed = True
			else:
				proceed = False

	elif command == 'setup_nmr':
		if stage in ['running-opt', 'pre-nmr']:
			proceed = True
		elif stage in ['init', 'pre-opt']:
			print('Not ready for NMR calculations, need to run init, conf_search, and setup_opt first')
			proceed = False
		else:
			print('Already setup NMR calculations, do you want to continue ? This may overwrite/remove existing com/log files')
			answer = yes_or_no()
			if answer:
				proceed = True
			else:
				proceed = False

	elif command == 'process_nmr':
		if stage in ['running-nmr']:
			proceed = True
		elif stage in ['init', 'pre-opt', 'running-opt']:
			print('Not ready to process NMR calculations, need to run init, conf_search, setup_opt, and setup_nmr first')
			proceed = False
		elif stage in ['post']:
			print('Re-running NMR processing')
			proceed = True

	elif command == 'update':
		if stage in ['init', 'pre-opt']:
			print('Not ready to update, must be running gaussian calculations, need to run init, conf_search and setup_opt first')
			proceed = False
		else:
			proceed = True

	elif command == 'check_status':
		proceed = True


	return proceed, molecule












##

from file_creation import HPC_submission as HPCsub
import file_read.orca_read as orcaread
import file_creation.orca_submission as orcasub


def setup_optimisation(molecule, prefs, path='', max=50):
	opt_files = []
	for conformer in molecule.conformers:
		conformer.opt_in = orcasub.make_optin(prefs, conformer, path + 'optimisation/')
		conformer.opt_log = conformer.optin.split('.')[0] + '.log'
		conformer.opt_status = 'pre-submission'
		opt_files.append(conformer.opt_in)

	IN_ARRAY = 'optimisation/OPT_IN_ARRAY.txt'
	with open(path + IN_ARRAY, 'w') as f:
		for file in opt_files:
			print(file, file=f)

	header = HPCsub.make_HPC_header(prefs)

	files = len(opt_files)
	chunks = HPCsub.get_chunks(files)
	for ck in range(chunks):
		start = (ck * max) + 1
		end = ((ck + 1) * max)
		if end > files:
			end = files
		strings = HPCsub.make_HPC_orca_batch_submission(prefs, IN_ARRAY, start, end, ck)

		if prefs['comp']['system'] == 'BC3':
			filename = path + 'OPT_' + molecule.molid + '_' + str(ck) + '.qsub'
		elif prefs['comp']['system'] == 'BC4':
			filename = path + 'OPT_' + molecule.molid + '_' + str(ck) + '.slurm'
		elif prefs['comp']['system'] == 'localbox':
			filename = path + 'OPT_' + molecule.molid + '_' + str(ck) + '.sh'
		with open(filename, 'w') as f:
			for string in header:
				print(string, file=f)
			for string in strings:
				print(string, file=f)

	print('Created ', chunks, ' submission files. . .')
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


def process_optimisation(molecule, prefs, path=''):

	good = 0
	bad = 0

	for conformer in molecule.conformers:
		status = orcaread.get_opt_status(conformer.opt_log)
		if status == 'successful':
			good +=1
			conformer.read_opt(file, type='orca')
		else:
			bad += 1

		if conformer.opt_status != status:
			conformer.nmr_status = 'None'

		conformer.opt_status = status

		string = 'Conformer {molid:^3s} status: {status:^10s} | Energy {energy:<10.4f}'.format(molid=str(conformer.molid),
																								status=conformer.opt_status,
																								energy=conformer.energy)
		print(string)

	print(good, ' successful optimisations, ', bad, ' failed, out of ', len(statuss))



































###

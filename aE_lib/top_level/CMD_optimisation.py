import sys
from file_creation import HPC_submission as HPCsub
import file_read.orca_read as orcaread
import file_creation.orca_submission as orcasub


def setup_optimisation(molecule, prefs, path='', max=50):
	opt_files = []

	for conformer in molecule.conformers:
		conformer.opt_in = orcasub.make_optin(prefs, conformer, path + 'optimisation/')
		conformer.opt_log = conformer.opt_in.split('.')[0] + '.log'
		conformer.opt_status = 'pre-submission'
		opt_files.append('optimisation/' + conformer.opt_in.split('/')[-1])

	if len(opt_files) == 0:
		print('No files to submit. . .')
		sys.exit(0)
	IN_ARRAY = 'optimisation/OPT_IN_ARRAY.txt'
	with open(path + IN_ARRAY, 'w') as f:
		for file in opt_files:
			print(file, file=f)

	system = prefs['comp']['system']
	memory = prefs['optimisation']['memory']
	processors = prefs['optimisation']['processors']
	walltime = prefs['optimisation']['walltime']

	files = len(opt_files)
	chunks = HPCsub.get_chunks(files)
	qsub_names = []
	for ck in range(chunks):
		start = (ck * max) + 1
		end = ((ck + 1) * max)
		if end > files:
			end = files


		#header = HPCsub.make_HPC_header(jobname=jobname, system=system, nodes=1, ppn=processors, walltime=walltime, mem=memory)
		jobname = 'aE_' + molecule.molid + '_' + str(ck) + '_optimisation'
		strings = HPCsub.make_HPC_orca_batch_submission(prefs, molecule.molid, IN_ARRAY, start, end, ck=ck, jobname=jobname, nodes=1, ppn=processors, walltime=walltime, mem=memory)

		if prefs['comp']['system'] == 'BC3':
			filename = path + 'OPT_' + molecule.molid + '_' + str(ck) + '.qsub'
		elif prefs['comp']['system'] == 'BC4':
			filename = path + 'OPT_' + molecule.molid + '_' + str(ck) + '.slurm'
		elif prefs['comp']['system'] == 'localbox':
			filename = path + 'OPT_' + molecule.molid + '_' + str(ck) + '.sh'
		with open(filename, 'w') as f:
			for string in strings:
				print(string, file=f)
		qsub_names.append(filename)

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

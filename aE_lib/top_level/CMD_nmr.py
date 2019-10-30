from file_creation import HPC_submission as HPCsub
import file_read.orca_read as orcaread
import file_creation.orca_submission as orcasub
import file_creation.structure_formats.nmredata as nmredata

def setup_nmr(molecule, prefs, path='', ids=[]):

	nmr_files = []
	for conformer in molecule.conformers:
		if conformer.nmr_status == 'None' and conformer.opt_status == 'successful':
			conformer.nmr_in = orcasub.make_nmrin(prefs, conformer, path + 'nmr/')
			conformer.nmr_log = conformer.nmr_in.split('.') + '.log'
			conformer.nmr_status = 'pre-submission'
			nmr_files.append(conformer.nmr_in)

	IN_ARRAY = 'nmr/NMR_IN_ARRAY.txt'
	with open(path + IN_ARRAY, 'w') as f:
		for file in opt_files:
			print(file, file=f)

	header = HPCsub.make_HPC_header(prefs)

	files = len(nmr_files)
	chunks = HPCsub.get_chunks(files)
	for ck in chunks:
		start = (ck * max) + 1
		end = ((ck + 1) * max)
		if end > files:
			end = files
		strings = HPCsub.make_orca_batch_submission(prefs, IN_ARRAY, start, end, ck)

		if prefs['comp']['system'] == 'BC3':
			filename = path + 'NMR_' + molecule.molid + '_' + str(ck) + '.qsub'
		elif prefs['comp']['system'] == 'BC4':
			filename = path + 'NMR_' + molecule.molid + '_' + str(ck) + '.slurm'
		elif prefs['comp']['system'] == 'localbox':
			filename = path + 'NMR_' + molecule.molid + '_' + str(ck) + '.sh'
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


def process_nmr(molecule, prefs, path=''):

	good = 0
	bad = 0
	process = False

	for conformer in molecule.conformers:
		status = orcaread.get_nmr_status(conformer.nmr_log)
		if status == 'successful':
			good +=1
			conformer.read_nmr(file, type='orca')
		else:
			bad += 1

		if conformer.nmr_status != status:
			process = True
			outname = path + 'OUTPUT/' + molecule.molid + '_' + conformer.molid + '.nmredata.sdf'
			nmredata.nmrmol_to_nmredata(molecule, outname)

		conformer.nmr_status = status
		string = 'Conformer {molid:^3s} status: {status:^10s}'.format(molid=str(conformer.molid),
																								status=conformer.nmr_status)
		print(string)

	if process:
		molecule.boltzmann_average()
		outname = path + 'OUTPUT/' + molecule.molid + '.nmredata.sdf'
		nmredata.nmrmol_to_nmredata(molecule, outname)

	print(good, ' successful NMR calculations, ', bad, ' failed, out of ', len(statuss))















###

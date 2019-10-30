import file_read.orca_read as orcaread
import file_creation.orca_submission as orcasub
from file_creation import HPC_submission as HPCsub

def setup_resubmission(molecule, prefs, path=args.path):

	for conformer in molecule.conformers:
		optstatus = orcaread.get_opt_status(conformer.opt_log)

		if optstatus != 'successful':
			conformer.opt_in = orcasub.make_optin(prefs, conformer, path + 'optimisation/')
			conformer.opt_log = conformer.optin.split('.')[0] + '.log'
			conformer.opt_status = 'pre-submission'
			opt_files.append(conformer.opt_in)

		nmrstatus = orcaread.get_nmr_status(conformer.nmr_log)
		if conformer.nmr_status != 'successful' and conformer.opt_status == 'successful':
			conformer.nmr_in = orcasub.make_nmrin(prefs, conformer, path + 'nmr/')
			conformer.nmr_log = conformer.nmr_in.split('.') + '.log'
			conformer.nmr_status = 'pre-submission'
			nmr_files.append(conformer.nmr_in)

	for tag, in_files in zip(['OPT', 'NMR'], [opt_files, nmr_files])

		files = len(in_files)
		chunks = HPCsub.get_chunks(files)
		for ck in chunks:
			start = (ck * max) + 1
			end = ((ck + 1) * max)
			if end > files:
				end = files
			strings = HPCsub.make_orca_batch_submission(prefs, in_files, start, end, ck)

			if prefs['comp']['system'] == 'BC3':
				filename = path + 'RESUB_' + tag + '_' + molecule.molid + '_' + str(ck) + '.qsub'
			elif prefs['comp']['system'] == 'BC4':
				filename = path + 'RESUB_' + tag + '_' molecule.molid + '_' + str(ck) + '.slurm'
			elif prefs['comp']['system'] == 'localbox':
				filename = path + 'RESUB_' + tag + '_' molecule.molid + '_' + str(ck) + '.sh'
			with open(filename, 'w') as f:
				for string in header:
					print(string, file=f)
				for string in strings:
					print(string, file=f)

		print('Created ', len(chunks), ' ',  tag 'resubmission files. . .')


















###

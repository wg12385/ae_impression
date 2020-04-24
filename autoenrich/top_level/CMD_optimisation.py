# Copyright 2020 Will Gerrard
#This file is part of autoenrich.

#autoenrich is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoenrich is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.

import sys
from autoenrich.file_creation import HPC_submission as HPCsub
import autoenrich.file_read.orca_read as orcaread
import autoenrich.file_read.g09_read as g09read
import autoenrich.file_read.g16_read as g16read
import autoenrich.file_creation.orca_submission as orcasub
import autoenrich.file_creation.g09_submission as g09sub
import autoenrich.file_creation.g16_submission as g16sub
import os
import glob


def setup_optimisation(molecule, prefs, path='', max=50):
	opt_files = []

	for conformer in molecule.conformers:

		if conformer.redundant:
			continue

		if prefs['optimisation']['software'] == 'orca':
			conformer.opt_in = orcasub.make_optin(prefs, conformer.molid, conformer.xyz, conformer.types,
														path + 'optimisation/')
		elif prefs['optimisation']['software'] == 'g09':
			conformer.opt_in = g09sub.make_optcom(prefs, conformer.molid, conformer.xyz, conformer.types,
														path + 'optimisation/')
		elif prefs['optimisation']['software'] == 'g16':
			conformer.opt_in = g16sub.make_optcom(prefs, conformer.molid, conformer.xyz, conformer.types,
														path + 'optimisation/')


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
	software = prefs['optimisation']['software']

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
		strings = HPCsub.make_HPC_batch_submission(prefs, molecule.molid, IN_ARRAY, start, end, software=software,
									jobname=jobname, nodes=1, ppn=processors, walltime=walltime, mem=memory)

		if prefs['comp']['system'] == 'PBS':
			filename = path + 'OPT_' + molecule.molid + '_' + str(ck) + '.qsub'
		elif prefs['comp']['system'] == 'slurm':
			filename = path + 'OPT_' + molecule.molid + '_' + str(ck) + '.slurm'
		elif prefs['comp']['system'] == 'local':
			filename = path + 'OPT_' + molecule.molid + '_' + str(ck) + '.sh'
		with open(filename, 'w') as f:
			for string in strings:
				print(string, file=f)
		qsub_names.append(filename)

	print('Created ', chunks, ' submission files. . .')
	if prefs['comp']['system'] == 'PBS':
		print('Submit the calculations using:')
		for file in qsub_names:
			print('qsub ', file)
	elif prefs['comp']['system'] == 'PBS':
		print('Submit the calculations using:')
		for file in qsub_names:
			print('bash ', file)
	elif prefs['comp']['system'] == 'slurm':
		print('Havent finished this yet, good luck pal. . . .')


def process_optimisation(molecule, prefs, path=''):

	good = 0
	bad = 0

	for conformer in molecule.conformers:
		if prefs['optimisation']['software'] == 'orca':
			status = orcaread.get_opt_status(conformer.opt_log)
		elif prefs['optimisation']['software'] == 'g09':
			status = g09read.get_opt_status(conformer.opt_log)
		elif prefs['optimisation']['software'] == 'g16':
			status = g16read.get_opt_status(conformer.opt_log)

		if status == 'successful':
			good +=1
			conformer.read_structure(conformer.opt_log.split('.')[0] + '.xyz', type='xyz')
			conformer.read_opt(conformer.opt_log, type=prefs['optimisation']['software'])
			#files_to_delete = glob.glob(conformer.opt_log.split('.')[0] + '*.tmp')
			#for file in files_to_delete:
			#	os.remove(file)

		else:
			bad += 1

		if conformer.opt_status != status:
			conformer.nmr_status = 'None'

		conformer.opt_status = status

		string = 'Conformer {molid:^3s} status: {status:^10s} | Energy {energy:<10.4f}'.format(molid=str(conformer.molid),
																								status=conformer.opt_status,
																								energy=conformer.energy)
		print(string)

	print(good, ' successful optimisations, ', bad, ' failed, out of ', good+bad)



































###

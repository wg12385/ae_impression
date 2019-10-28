# functions for creating gaussian com files
from reference.periodic_table import Get_periodic_table



def make_optin(prefs, molname, xyz, type, directory=''):

	charge = prefs['mol']['charge']
	multiplicity = prefs['mol']['multiplicity']
	functional = prefs['optimisation']['functional']
	basis_set = prefs['optimisation']['basisset']
	solvent = prefs['optimisation']['solvent']
	direct_cmd_line_opt = prefs['optimisation']['custom_cmd_line']

	Periodic_table = Get_periodic_table()

	instr = '! ' + str(functional) + ' ' + str(basis_set) + ' TightSCF OPT'

	if solvent != 'none':
		instr += ' CPCM(' + solvent + ')'

	if direct_cmd_line_opt:
		instr1 = direct_cmd_line_opt

	infile = directory.strip() + molname.strip() + '_OPT.in'

	strings = []
	strings.append(instr)
	strings.append('')
	strings.append("* xyz {0:<1d} {1:<1d}".format(charge, multiplicity))

	for i in range(len(xyz)):
		str_type = Periodic_table[type[i]]
		string = " {0:<2s}        {1:>10.5f}        {2:>10.5f}        {3:>10.5f}".format(str_type, xyz[i][0], xyz[i][1], xyz[i][2])
		strings.append(string)

	strings.append('*')
	strings.append('')
	strings.append('%geom')
	strings.append('     AddExtraBonds true         # switch on/off assigning bonds to atom pairs that are')
	strings.append('                                #  connected by more than <Max_Length> bonds and are less')
	strings.append('                                #  than <MaxDist> Ang. apart (default true)')
	strings.append('     AddExtraBonds_MaxLength 10 # cutoff for number of bonds connecting the two')
	strings.append('                                #  atoms (default 10)')
	strings.append('     AddExtraBonds_MaxDist 5    # cutoff for distance between two atoms (default 5 Ang.)')
	strings.append('end')

	with open(infile, 'w') as f_handle:
		for string in strings:
			print(string, file=f_handle)

	return infile

def make_nmrin(prefs, molname, xyz, type, directory=''):

	charge = prefs['mol']['charge']
	multiplicity = prefs['mol']['multiplicity']
	functional = prefs['NMR']['functional']
	basis_set = prefs['NMR']['basisset']
	aux_basis_set = prefs['NMR']['aux_basis_set']
	solvent = prefs['NMR']['solvent']
	direct_cmd_line_opt = prefs['NMR']['custom_cmd_line']


	Periodic_table = Get_periodic_table()
	try:
		int(memory)
		int(processors)
	except:
		return

	instr = '! ' + str(functional) + ' ' + str(basis_set) + ' ' + str(aux_basis_set) + 'TightSCF' + ' NMR '

	if solvent != 'none':
		instr += ' CPCM(' + solvent + ')'

	if direct_cmd_line_nmr:
		instr = direct_cmd_line_nmr

	infile = directory.strip() + molname.strip() + '_NMR.in'

	strings = []
	strings.append(instr)
	strings.append("")
	strings.append("* xyz {0:<1d} {1:<1d}".format(charge, multiplicity))
	for i in range(len(xyz)):
		str_type = Periodic_table[type[i]]
		string = " {0:<2s}        {1:>10.6f}        {2:>10.6f}        {3:>10.6f}".format(str_type, xyz[i][0], xyz[i][1], xyz[i][2])
		strings.append(string)
	strings.append('*')

	strings.append('%epnmr')
	for type in prefs['NMR']['shift_nuclei']:
		strings.append('       Nuclei = all {type:<2s} { shift }'.format(Periodic_table(type)))
	for type in prefs['NMR']['spin_nuclei']:
		strings.append('       Nuclei = all {type:<2s} {ssall}'.format(Periodic_table(type)))
	strings.append('SpinSpinRThresh {0:<f}'.format(prefs['NMR']['spin_thresh']))
	strings.append('end')

	with open(infile, 'w') as f_handle:
		for string in strings:
			print(string, file=f_handle)

	return infile

def make_submission_array(molname, innames, path='', tag=''):
	filename = path + molname + '_' + tag + '_in_ARRAY'
	with open(filename, 'w') as f:
		for name in innames:
			name = name.split('/')[-1]
			print(name, file=f)

	return filename

def make_submission_qsub(prefs, in_array, innames, molname, start=-1, end=-1, path='', tag='',
							nodes=1, ppn=8, walltime='100:00:00', mem=32,):

	system = prefs['comp']['system']
	parallel = prefs['comp']['parallel']


	max = 50
	files = len(innames)
	if start < 0 and end < 0:
		start = 1
		if files > max:
			chunks = files / max
			if files % max > 0:
				chunks = int(chunks) + 1
			if chunks == 0:
				chunks = 1
		else:
			chunks = 1
			end = files
	else:
		if end - start >= max:
			chunks = files / max
			if files % max > 0:
				chunks = int(chunks) + 1
		else:
			chunks = 1

	chunks = int(chunks)
	filenames = []
	for ck in range(chunks):
		if system == 'BC3':
			filename = path + molname + '_' + tag + '_' + str(ck) + '.qsub'
		elif system == 'Grendel':
			filename = path + 'submit_' + molname + '_' + tag + '_' + str(ck) + '.sh'
		elif system == 'BC4':
			filename = path + molname + '_' + tag + '_' + str(ck) + '.slurm'
		elif system == 'localbox':
			filename = path + molname + '_' + tag + '_' + str(ck) + '.sh'
		filenames.append(filename)
		start = (ck * max) + 1
		end = ((ck + 1) * max)
		if end > files:
			end = files
		with open(filename, 'w') as f:
			strings = []
			strings.append('#!/bin/bash')
			if system == 'BC3':
				strings.append("#PBS -l nodes={0:<1d}:ppn={1:<1d}".format(nodes, ppn))
				strings.append("#PBS -l walltime={0:<9s}".format(walltime))
				strings.append("#PBS -l mem={0:<1d}GB".format(mem))
				if parallel:
					strings.append("#PBS -N {0:>1s}_{1:<1d}".format(molname, ck))
					strings.append("#PBS -t {0:>1d}-{1:<1d}".format(start, end))
					strings.append("cd $PBS_O_WORKDIR")
					strings.append("NMRNAME=$(gawk -v y=${{PBS_ARRAYID}} 'NR == y' {0:<5s})".format(in_array))
					strings.append("g09 $NMRNAME")
				else:
					strings.append("cd $PBS_O_WORKDIR")
					strings.append("for i in $(seq {0:>1d} 1 {1:<1d});".format(start, end))
					strings.append("do")
					strings.append("NMRNAME=$(gawk -v y=${i} 'NR == y' {0:<5s})".format(in_array))
					strings.append("base=${NMRNAME::${#NMRNAME}-4}")
					strings.append("FILE=${base}.log")
					strings.append("if test -f '$FILE'; then")
					strings.append("continue")
					strings.append("else")
					strings.append("g09 $NMRNAME")
					strings.append("fi")
					strings.append("done")

			elif system == 'Grendel':
				strings.append("NUMBERS=$(seq {0:>1d} {1:<1d})".format(start, end))
				strings.append("for NUM in ${NUMBERS}; do")
				strings.append("  NMRNAME=$(head -n${{NUM}} {0:<5s} | tail -1)".format(in_array))
				strings.append("  orca  $NMRNAME -l nodes={0:<1d}:ppn={1:<1d} -l walltime={2:<9s} -l mem={3:<1d}GB -N {4:>1s}_{5:<1d}_${{NUM}}".format(nodes, ppn, walltime, mem, molname, ck))
				strings.append("  done")
			elif system == 'BC4':
				print('not done yet . . .')
			elif system == 'localbox':
				strings.append("NUMBERS=$(seq {0:>1d} {1:<1d})".format(start, end))
				strings.append("for NUM in ${NUMBERS}; do")
				strings.append("  NMRNAME=$(head -n${{NUM}} {0:<5s} | tail -1)".format(com_array))
				strings.append("  OUTNAME=$( echo $NMRNAME | sed 's/.com/.log/')")
				strings.append("  orca ${NMRNAME} > ${OUTNAME}")
				strings.append("done")
			else:
				print('preference value for system: ', system, ' was not recognised, accepted values are BC3, BC4, Grendel, localbox')

			for string in strings:
				print(string, file=f)

	return filenames

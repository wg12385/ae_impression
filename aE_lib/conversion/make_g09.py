# functions for creating gaussian com files
from reference.periodic_table import Get_periodic_table



def make_optcom(molname, xyz, type, directory='',
					charge=0, multiplicity=1, memory=26, processors=8,
					opt='tight', freq=True, functional='mpw1pw91', basis_set='6-31g(d)',
					solvent='none', solventmodel='iefpcm' , grid='fine', direct_cmd_line_opt = 'none'):
	## Inputs:
	#  molname      = molecule name, string. 'Progesterone'
	#  xyz          = xyz coordinates, Nx3 numpy array.
	#  type         = atom types, list(or 1D array) length N.
	#  charge       = molecule charge, integer.
	#  multiplicity = multiplcity of molecule, integer.
	#  memory       = memory value to use for gaussian, integer.
	#  processors   = number of processoers to use for gaussian, integer.
	#  instr        = instruction line for gaussian, string. 'opt=tight freq mpw1pw91/6-311g(d,p) scrf=(iefpcm,solvent=chloroform) geom=distance MaxDisk=50GB'
	Periodic_table = Get_periodic_table()
	try:
		int(memory)
		int(processors)
	except:
		return

	if freq == "True":
		instr = 'opt=' + str(opt) + ' freq ' + str(functional) + '/' + str(basis_set) + ' integral=' + str(grid) + ' MaxDisk=50GB'
	else:
		instr = 'opt=' + str(opt) + ' ' + str(functional) + '/' + str(basis_set) + ' integral=' + str(grid) + ' MaxDisk=50GB'

	if solvent != 'none':
		instr += ' scrf=(' + str(solventmodel) + ',solvent=' + str(solvent) + ')'

	if direct_cmd_line_opt != 'none':
		instr1 = direct_cmd_line_opt

	comfile = directory.strip() + molname.strip() + '_OPT.com'


	with open(comfile, 'w') as f_handle:
		strings = []
		strings.append("%mem={0:<1d}GB".format(memory))
		strings.append("%NProcShared={0:<1d}".format(processors))
		strings.append("# {0:<1s}".format(instr))
		strings.append("")
		strings.append("{0:<1s} OPT".format(molname))
		strings.append("")
		strings.append("{0:<1d} {1:<1d}".format(charge, multiplicity))
		for string in strings:
			print(string, file=f_handle)
		for i in range(len(xyz)):
			str_type = Periodic_table[type[i]]
			string = " {0:<2s}        {1:>10.5f}        {2:>10.5f}        {3:>10.5f}".format(str_type, xyz[i][0], xyz[i][1], xyz[i][2])
			print(string, file=f_handle)
		for i in range(4):
			print("", file=f_handle)
		strings = []

	return comfile

def make_nmrcom(molname, xyz, type, directory='',
				charge=0, multiplicity=1, memory=26, processors=8,
				functional='wb97xd', basis_set='6-311g(d,p)', mixed=True,
				solvent='none', solventmodel='iefpcm',  direct_cmd_line_nmr = 'none'):

	## Inputs:
	#  molname      = molecule name, string. 'Progesterone'
	#  xyz          = xyz coordinates, Nx3 numpy array.
	#  type         = atom types, list(or 1D array) length N.
	#  dir          = directory where optimisation chk files are stored, string. 'progesterone_opt_chk_files/'
	#  charge       = molecule charge, integer.
	#  multiplicity = multiplcity of molecule, integer.
	#  memory       = memory value to use for gaussian, integer.
	#  processors   = number of processoers to use for gaussian, integer.
	#  instr        = instruction line for gaussian, string. 'opt=tight freq mpw1pw91/6-311g(d,p) scrf=(iefpcm,solvent=chloroform) geom=distance MaxDisk=50GB'
	Periodic_table = Get_periodic_table()
	try:
		int(memory)
		int(processors)
	except:
		return

	if mixed == "True":
		instr='nmr(giao,spinspin,mixed)' + str(functional) + '/' + str(basis_set) + ' maxdisk=50GB'
	else:
		instr='nmr(giao,spinspin)' + str(functional) + '/' + str(basis_set) + ' maxdisk=50GB'

	if solvent != 'none':
		instr += ' scrf=(' + str(solventmodel) + ',solvent=' + str(solvent) + ')'

	if direct_cmd_line_nmr != 'none':
		instr = direct_cmd_line_nmr

	comfile = directory.strip() + molname.strip() + '_NMR.com'
	chkfile = molname.strip() + '_NMR'
	with open(comfile, 'w') as f_handle:
		strings = []
		strings.append("%Chk={0:<1s}".format(chkfile))
		strings.append("%NoSave")
		strings.append("%mem={0:<1d}GB".format(memory))
		strings.append("%NProcShared={0:<1d}".format(processors))
		strings.append("#T {0:<1s}".format(instr))
		strings.append("")
		strings.append("{0:<1s} NMR".format(molname))
		strings.append("")
		strings.append("{0:<1d} {1:<1d}".format(charge, multiplicity))
		for string in strings:
			print(string, file=f_handle)
		for i in range(len(xyz)):
			str_type = Periodic_table[type[i]]
			string = " {0:<2s}        {1:>10.6f}        {2:>10.6f}        {3:>10.6f}".format(str_type, xyz[i][0], xyz[i][1], xyz[i][2])
			print(string, file=f_handle)
		for i in range(4):
			print("", file=f_handle)

	return 0

def make_submission_array(molname, comnames, path='', tag=''):
	filename = path + molname + '_' + tag + '_COM_ARRAY'
	with open(filename, 'w') as f:
		for name in comnames:
			print(name, file=f)

	return filename

def make_submission_qsub(com_array, comnames, molname, path='', parallel=False, system='BC3',
							tag='', nodes=1, ppn=8, walltime='100:00:00',
								mem=32, start=-1, end=-1):

	max = 50
	files = len(comnames)
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
		else:
			filename = path + 'submit_' + molname + '_' + tag + '_' + str(ck) + '.sh'
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
					strings.append("NMRNAME=$(gawk -v y=${{PBS_ARRAYID}} 'NR == y' {0:<5s})".format(com_array))
					strings.append("g09 $NMRNAME")
				else:
					strings.append("cd $PBS_O_WORKDIR")
					strings.append("for i in $(seq {0:>1d} 1 {1:<1d});".format(start, end))
					strings.append("do")
					strings.append("NMRNAME=$(gawk -v y=${i} 'NR == y' {0:<5s})".format(com_array))
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
				strings.append("  NMRNAME=$(head -n${{NUM}} {0:<5s} | tail -1)".format(com_array))
				strings.append("  qg09  $NMRNAME -l nodes={0:<1d}:ppn={1:<1d} -l walltime={2:<9s} -l mem={3:<1d}GB -N {4:>1s}_{5:<1d}_${{NUM}}".format(nodes, ppn, walltime, mem, molname, ck))
				strings.append("  done")
			else:
				print('preference value for system: ', system, ' was not recognised, accepted values are BC3 or Grendel')

			for string in strings:
				print(string, file=f)

	return filenames

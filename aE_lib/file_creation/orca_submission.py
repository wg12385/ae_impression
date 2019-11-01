# make orca submission files

# functions for creating gaussian com files
from reference.periodic_table import Get_periodic_table

def make_optin(prefs, mol, directory=''):

	molname = mol.molid
	xyz = mol.xyz
	types = mol.types

	charge = prefs['mol']['charge']
	multiplicity = prefs['mol']['multiplicity']
	functional = prefs['optimisation']['functional']
	basis_set = prefs['optimisation']['basisset']
	solvent = prefs['optimisation']['solvent']
	direct_cmd_line_opt = prefs['optimisation']['custom_cmd_line']

	Periodic_table = Get_periodic_table()

	instr = '! ' + str(functional) + ' ' + str(basis_set) + ' TightSCF OPT miniprint'

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
		str_type = Periodic_table[types[i]]
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

def make_nmrin(prefs, mol, directory=''):

	molname = mol.molid
	xyz = mol.xyz
	types = mol.types

	charge = prefs['mol']['charge']
	multiplicity = prefs['mol']['multiplicity']
	functional = prefs['NMR']['functional']
	basis_set = prefs['NMR']['basisset']
	aux_basis_set = prefs['NMR']['aux_basis_set']
	solvent = prefs['NMR']['solvent']
	direct_cmd_line_nmr = prefs['NMR']['custom_cmd_line']


	Periodic_table = Get_periodic_table()

	instr = '! ' + str(functional) + ' ' + str(basis_set) + ' ' + str(aux_basis_set) + 'TightSCF miniprint' + ' NMR '

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
		strings.append("       Nuclei = all {type:<2s}".format(type=Periodic_table[type]) + '  { shift }')
	for type in prefs['NMR']['spin_nuclei']:
		strings.append("       Nuclei = all {type:<2s}".format(type=Periodic_table[type]) + '  { ssall }')
	strings.append('SpinSpinRThresh {0:<f}'.format(prefs['NMR']['spin_thresh']))
	strings.append('end')

	with open(infile, 'w') as f_handle:
		for string in strings:
			print(string, file=f_handle)

	return infile

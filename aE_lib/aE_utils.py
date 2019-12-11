# master script to do everything IMPRESSION-y

import argparse
import sys
import json
sys.path.append('aE_lib/')

from util.flag_handler import hdl_targetflag, flag_combos
from util.argparse_wizard import run_wizard

from top_level import CMD_trainmodel, CMD_predict, CMD_compare

from util.header import print_header_IMP

from molecule.nmrmol import nmrmol
from file_creation.structure_formats import nmredata
from file_read.g09_read import read_functional
from reference.tantillo import Get_tantillo_factors
import glob


if __name__ == "__main__":
	# Argparser
	parser = argparse.ArgumentParser(description='aE_utils')

	parser.add_argument('Command', help='aE_utils command',
							choices=['convert_to_nmredata'])

	parser.add_argument('--files', help='Files to process',
							action="store", dest='files', default='')

	parser.add_argument('--type', help='File type',
							action="store", dest='type', default='g09')

	parser.add_argument('--out_dir', help='output directory',
							action="store", dest='out_dir', default='')

	args = vars(parser.parse_args())

	if args['Command'] == 'convert_to_nmredata':

		for f, file in enumerate(glob.glob(args['files'])):

			print(file)
			mol = nmrmol(molid=f)
			#mol.read_structure(file, args['type'])
			mol.read_nmr(file, args['type'])

			'''
			if args['type'] in ['orca', 'g09']:
				#functional, basis_set = read_functional(file)
				functional = 'wb97xd'
				basis_set = '6311gdp'
				scaling_factors = Get_tantillo_factors(basis_set, functional)
				mol.scale_shifts(scaling_factors)
			'''

			outname = file.split('/')[-1].split('.')[0]
			outfile = args['out_dir'] + outname + '.nmredata.sdf'
			nmredata.nmrmol_to_nmredata(mol, outfile)

















#

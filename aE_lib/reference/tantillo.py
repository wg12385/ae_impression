
import json
import re
import sys
import os

def Get_tantillo_factors(basis_set='6-311g(d,p)', functional='wb97xd'):
	'''
	tantillo=	{}
	tantillo['functional'] = functional
	tantillo['basis_set'] = basis_set
	tantillo['H'] = [-1.0719, 32.1254]
	tantillo['C'] = [-1.0399, 187.136]

	dir_path = os.path.dirname(os.path.realpath(__file__))
	json_file = dir_path + '/scaling_factors/' + re.sub(r'\W+', '', functional+basis_set) + '.json'
	json.dump(tantillo, open(json_file, 'w'), indent=4)
	'''
	try:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		json_file = dir_path + '/scaling_factors/' + re.sub(r'\W+', '', functional+basis_set) + '.json'
		tantillo = json.load(open(json_file, 'r'))
	except Exception as e:
		print(e)
		print('Couldnt load scaling factor json file. . . ')
		sys.exit(0)

	return tantillo



def Calculate_new_factors(exp_mols, dft_mols, functional, basis_set):
	print('STUB, havent written this yet')
	# do basic linear regression


	dir_path = os.path.dirname(os.path.realpath(__file__))
	json_file = dir_path + '/scaling_factors/' + re.sub(r'\W+', '', functional+basis_set) + '.json'
	json.dump(tantillo, open(json_file, 'w'), indent=4)






















#

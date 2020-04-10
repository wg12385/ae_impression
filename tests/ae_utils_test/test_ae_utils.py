import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../../')


from autoenrich.molecule.nmrmol import nmrmol
from autoenrich.top_level import CMD_compare

import glob


def test_convertnmredata():
	status = "SUCCESS"



	return status




def test_comparedatasets():
	status = 'SUCCESS'

	path = '/'.join(__file__.split('/')[:-1])

	case = {}
	case['comp_sets'] = [path + '/test_files/IMP*.sdf', path + '/test_files/DFT*.sdf']
	case['match_criteria'] = 'id'
	case['comp_targets'] = ['HCS', 'CCS', '1JCH', '3JHH']
	case['comp_labels'] = ['IMP', 'DFT']
	case['output_path'] = path + '/../test_tmp'

	CMD_compare.compare_datasets(case)

	files_to_remove = glob.glob(path + '/../test_tmp/*')
	for file in files_to_remove:
		os.remove(file)


	return status




if __name__ == "__main__":

	status1 = test_convertnmredata()
	status2 = test_comparedatasets()


	print(status1, status2)

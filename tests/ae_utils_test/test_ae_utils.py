from autoENRICH.molecule.nmrmol import nmrmol
from autoENRICH.top_level import CMD_compare
print(CMD_compare.__file__)


def test_convertnmredata():
	status = "SUCCESS"



	return status




def test_comparedatasets():
	status = 'SUCCESS'

	case = {}
	case['comp_sets'] = ['test_files/IMP*.sdf', 'test_files/DFT*.sdf']
	case['match_criteria'] = 'id'
	case['comp_targets'] = ['HCS', 'CCS', '1JCH', '3JHH']
	case['comp_labels'] = ['IMP', 'DFT']

	CMD_compare.compare_datasets(case)


	return status




if __name__ == "__main__":

	status1 = test_convertnmredata()
	status2 = test_comparedatasets()


	print(status1, status2)

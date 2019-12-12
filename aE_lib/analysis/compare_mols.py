
from util.flag_handler.hdl_targetflag import target_to_flag
import numpy as np

def mol_isequal(mol1, mol2, print=False):

	if len(mol1.types) != len(mol2.types):
		if print:
			print('Molecules are different sizes')
		return False
	if len(mol1.xyz) != len(mol2.xyz):
		if print:
			print('Molecules are different sizes')
		return False

	atoms = len(mol1.types)

	# Check for same chemical formula
	counts1 = np.zeros(1000, dtype=np.int32)
	counts2 = np.zeros(1000, dtype=np.int32)
	for i in range(atoms):
		counts1[mol1.types[i]] += 1
		counts2[mol2.types[i]] += 1
	Fail = False
	for i in range(1000):
		if counts1[i] == 0 and counts2[i] == 0:
			continue
		if counts1[i] != counts2[i]:
			Fail = True
	if Fail:
		if print:
			print('Molecules have different compositions')
		return False

	# Check for same connectivity
	cpl_dict1 = {}
	cpl_dict2 = {}

	for i in range(atoms):
		for j in range(atoms):
			if mol1.conn[i][j] != 0:
				flag1 = target_to_flag([mol1.coupling_len[i][j], mol1.types[i], mol1.types[j]])
				flag2 = target_to_flag([mol2.coupling_len[i][j], mol2.types[i], mol2.types[j]])

				if flag1 not in cpl_dict1.keys():
					cpl_dict1[flag1] = 0
				if flag2 not in cpl_dict2.keys():
					cpl_dict2[flag2] = 0

				cpl_dict1[flag1] = cpl_dict1[flag1] + 1
				cpl_dict2[flag2] = cpl_dict2[flag2] + 1

	if cpl_dict1 != cpl_dict2:
		if print:
			print('Connectivity is different')
			return False

	return True


















###

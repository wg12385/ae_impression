from molecule.dataset import dataset
from util.flag_handler.hdl_targetflag import flag_to_target
from util.filename_utils import get_unique_part
from analysis.compare_mols import mol_isequal
from file_creation.data_formats.csv import print_mol_csv
import glob
import numpy as np

def compare_datasets(args):

	att_mols = []
	sets = []
	for set_list in args['comp_sets']:
		print('Getting molecules from ', set_list)
		set = dataset()
		label_part = get_unique_part(glob.glob(set_list))
		set.get_mols(glob.glob(set_list), label_part=label_part)
		print(len(set.mols), ' molecules found from ', len(glob.glob(set_list)), ' files')

		sets.append(set)

	assert len(sets) > 1, print('Only one set found. . .')
	assert len(sets[0].mols) == len(sets[1].mols), print('Different numbers of molecules in sets')

	for m1, mol1 in enumerate(sets[0].mols):
		for m2, mol2 in enumerate(sets[1].mols):

			if not mol_isequal(mol1, mol2):
				continue

			if [mol1, mol2] in att_mols:
				continue

			if len(sets) > 2:
				for m3, mol3 in enumerate(sets[2].mols):
					if not mol_isequal(mol1, mol3):
						continue

					if [mol1, mol2, mol3] in att_mols:
						continue

					att_mols.append([mol1, mol2, mol3])

			else:
				att_mols.append([mol1, mol2])

	print(len(att_mols), ' molecules matched, out of ', len(sets[0].mols))



	for targetflag in args['comp_targets']:
		print(targetflag)
		target = flag_to_target(targetflag)

		for set in sets:
			set.get_features_frommols({'featureflag': 'dummy',
										'targetflag': targetflag,
										'max_size': 0})

		values = []
		refs = []
		typerefs = []

		assert len(sets[0].r) == len(sets[1].r)
		if len(sets) > 2:
			assert len(sets[2].r) == len(sets[1].r)

		for group in att_mols:
			for i in range(len(sets[0].r)):
				ref1 = sets[0].r[i]
				if ref1[0] != group[0].molid:
					continue
				val1 = sets[0].y[i]
				typeref1 = [group[0].types[row-1] for row in ref1[1:]]

				for j in range(len(sets[0].r)):

					ref2 = sets[1].r[j]

					if ref2[0] != group[1].molid:
						continue
					val2 = sets[1].y[j]
					typeref2 = [group[1].types[row-1] for row in ref2[1:]]

					bad = False
					for xx in range(1, len(ref1)):
						if ref1[xx] != ref2[xx]:
							bad = True
					if bad:
						continue

					if len(sets) > 2:
						for k in range(len(sets[0].r)):
							ref3 = sets[2].r[k]
							if ref3[0] != group[2].molid:
								continue
							val3 = sets[2].y[k]
							typeref3 = [group[2].types[row-1] for row in ref3[1:]]


							bad = False
							for xx in range(1, len(ref1)):
								if ref1[xx] != ref3[xx]:
									bad = True
							if bad:
								continue

							refs.append([ref1, ref2, ref3])
							values.append([val1, val2, val3])
							typerefs.append([typeref1, typeref2, typeref3])
					else:
						refs.append([ref1, ref2])
						values.append([val1, val2])
						typerefs.append([typeref1, typeref2])

		outname = 'Comparison_' + str(targetflag) + '.csv'
		print_mol_csv(outname, refs, typerefs, values, args['comp_labels'])

		x = [row[0] for row in values]
		y = [row[1] for row in values]

		MAE = np.mean(np.absolute(np.asarray(x)-np.asarray(y)))
		print('MAE between ', args['comp_labels'][0], args['comp_labels'][1], ' = ', MAE)
















##

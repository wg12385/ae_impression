from ml.models import FCHLmodel, KRRmodel, TFMmodel
import numpy as np
import itertools
import pickle
import sys


def grid_search(x, y, modelflag='KRR', featureflag='CMAT', id='test_model', logfile='GS.log',
 				param_ranges={}, param_logs={}, grid_density=11, cv_steps=5):

	# create model
	if modelflag == 'KRR':
		model = KRRmodel.KRRmodel(id, x, y, params={})
	elif modelflag == 'FCHL':
		model = FCHLmodel.FCHLmodel(id, x, y, params={})
	elif modelflag == 'TFM':
		model = TFMmodel.TFMmodel(id, x, y, params={})

	# define grid
	search_grid = {}
	search_list = []
	size_list = []

	# determine whether log dictionary was provided
	if len(param_logs) == 0:
		check_logs = False
	else:
		check_logs = True

	# create search grid
	for param in param_ranges:
		if check_logs:
			if param_logs[param]:
				search_grid[param] = np.logspace(param_ranges[param][0], param_ranges[param][1], grid_density, endpoint=True)
			else:
				search_grid[param] = np.linspace(param_ranges[param][0], param_ranges[param][1], grid_density, endpoint=True)
		else:
			search_grid[param] = np.linspace(param_ranges[param][0], param_ranges[param][1], grid_density, endpoint=True)

	# create list of parameter sets to test
	for i in range(grid_density ** num_params):
		param_point = {}
		xx = 0
		for param in search_grid.keys():
			idx = int(i / grid_density ** xx) % grid_density
			param_point[param] = search_grid[param][idx]
			xx += 1
		search_list.append(param_point)


	strings = []
	strings.append('HPS GRID SEARCH')
	strings.append(modelflag + '   ' + featureflag + '   ' + paramflag)
	for param in param_ranges.keys():
		strings.append('{param:<10s}: {low:>10.4g}  <--->  {high:<10.4g}'.format(param=param,
																	low=param_ranges[param][0],
																	high=param_ranges[param][1]))
	strings.append('')
	strings.append('START')
	string = '{i:<10s}\t{score:<10s}'.format(i='i', score='SCORE')
	for param in param_ranges.keys():
		string = string + '\t{param:<15s}'.format(param)
	strings.append(string)
	with open(logfile, 'w') as f:
		for string in strings:
			print(string, file=f)

	BEST_SCORE = 999.9999999
	BEST_PARAMS = {}
	for p, params in enumerate(search_list):

		model.params = params

		y_pred = model.cv_predict(cv_steps)

		score = np.mean(np.absolute(y_pred - y))

		with open(logfile, 'a') as f:
			string = '{i:<10d}\t{score:<10.5f}'.format(i=p, score=score)
			for param in params.keys():
				string = string + '\t{param:<15.4g}'.format(params[param])

		if score < BEST_SCORE:
			BEST_SCORE = score
			BEST_PARAMS = params


	# create optimised model and save
	model.params = BEST_PARAMS
	model.train()

	outname = id + '_model.pkl'
	pickle.dump(open(outname, "wb"))

	print('Optimised model saved in ', outname)

	return score















###

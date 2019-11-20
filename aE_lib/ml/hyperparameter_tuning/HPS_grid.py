from ml.hyperparameter_tuning import HPS_generic as generic
import numpy as np
import itertools
import pickle
import sys

def grid_search(dataset, args):

	# define grid
	search_grid = {}
	search_list = []
	size_list = []

	# determine whether log dictionary was provided
	if len(args['param_logs']) == 0:
		check_logs = False
	else:
		check_logs = True

	# create search grid
	for param in args['param_ranges']:
		if check_logs:
			if args['param_logs'][param] == 'log':
				search_grid[param] = np.logspace(args['param_ranges'][param][0], args['param_ranges'][param][1], args['grid_density'], endpoint=True)
			else:
				search_grid[param] = np.linspace(args['param_ranges'][param][0], args['param_ranges'][param][1], args['grid_density'], endpoint=True)
		else:
			search_grid[param] = np.linspace(args['param_ranges'][param][0], args['param_ranges'][param][1], args['grid_density'], endpoint=True)

	# create list of parameter sets to test
	for i in range(args['grid_density'] ** len(args['param_ranges'].keys())):
		param_point = {}
		xx = 0
		for param in search_grid.keys():
			idx = int(i / args['grid_density'] ** xx) % args['grid_density']
			param_point[param] = search_grid[param][idx]
			xx += 1
		search_list.append(param_point)

	generic.setup_logfile(args)

	BEST_SCORE = 999.9999999
	BEST_PARAMS = {}
	for p, params in enumerate(search_list):

		BEST_SCORE, BEST_PARAMS = generic.HPS_iteration(p, dataset, args, next_point_to_probe=next_point_to_probe,
															BEST_SCORE=BEST_SCORE, BEST_PARAMS=BEST_PARAMS)

	generic.save_models(dataset, BEST_PARAMS, args)
	print('Optimised model saved in ', outname)

	return dataset, BEST_SCORE















###

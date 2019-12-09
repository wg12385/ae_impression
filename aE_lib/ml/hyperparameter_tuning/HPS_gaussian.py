
import bayes_opt as bayes
from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction
from bayes_opt.observer import JSONLogger
from bayes_opt.event import Events
from bayes_opt.util import load_logs

import pickle
import numpy as np

from ml.hyperparameter_tuning import HPS_generic as generic

def gaussian_search(dataset, args):

	# determine whether log dictionary was provided
	if len(args['param_logs']) == 0:
		check_logs = False
	else:
		check_logs = True

	pbounds = {}
	for param in args['param_ranges'].keys():
		pbounds[param] = (args['param_ranges'][param][0], args['param_ranges'][param][1])

	optimizer = BayesianOptimization(
		f=None,
		pbounds=pbounds,
		verbose=0, # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
		random_state=None
	)
	utility = UtilityFunction(kind="ucb", kappa=args['kappa'], xi=args['xi'])

	generic.setup_logfile(args)

	BEST_SCORE = 999.999
	BEST_PARAMS = {}
	searched = []
	for e in range(int(args['epochs'])):

		if args['random'] > 0 and e%args['random'] == 0:
			next_point_to_probe = {}
			for param in args['param_ranges'].keys():
				next_point_to_probe[param] = np.random.uniform(args['param_ranges'][param][0], args['param_ranges'][param][1])
		else:
			next_point_to_probe = optimizer.suggest(utility)

		not_unique = False
		attempts = 0
		while not_unique == False:
			switch = False
			attempts += 1
			for parms in searched:
				if next_point_to_probe == parms:
					switch = True
					next_point_to_probe = {}
					for param in args['param_ranges'].keys():
						next_point_to_probe[param] = np.random.uniform(args['param_ranges'][param][0], args['param_ranges'][param][1])
			if switch:
				not_unique == False
			elif attempts > 5:
				print('Search space expended for current parameters, finishing. . . ')
				outname = generic.save_models(dataset, BEST_PARAMS, args)
				print('Optimised model(s) saved in ', outname)
				return dataset, BEST_SCORE


		if check_logs:
			for param in args['param_ranges'].keys():
				if 'log' in args['param_logs'][param]:
					next_point_to_probe[param] = 10**next_point_to_probe[param]

		score, BEST_SCORE, BEST_PARAMS = generic.HPS_iteration(e, dataset, args, next_point_to_probe=next_point_to_probe,
															BEST_SCORE=BEST_SCORE, BEST_PARAMS=BEST_PARAMS)


		searched.append(next_point_to_probe)
		optimizer.register(params=next_point_to_probe, target=-score)



	outname = generic.save_models(dataset, BEST_PARAMS, args)
	print('Optimised model(s) saved in ', outname)

	return dataset, BEST_SCORE

























###

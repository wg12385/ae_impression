from ml.hyperparameter_tuning import HPS_generic as generic
from ml.models import FCHLmodel, KRRmodel, TFMmodel
import numpy as np
import itertools
import pickle
import sys


def random_search(dataset, args):

	# determine whether log dictionary was provided
	if len(args['param_logs']) == 0:
		check_logs = False
	else:
		check_logs = True


	generic.setup_logfile(args)

	BEST_SCORE = 999.999
	BEST_PARAMS = {}
	next_point_to_probe = {}

	for _ in range(args['epochs']):
		for param in args['param_ranges'].keys():
			next_point_to_probe[param] = np.random.uniform(args['param_ranges'][param][0], args['param_ranges'][param][1])

		if check_logs:
			for param in args['param_ranges'].keys():
				if args['param_logs'][param] == 'log':
					next_point_to_probe[param] = 10**next_point_to_probe[param]


		BEST_SCORE, BEST_PARAMS = generic.HPS_iteration(dataset, args, params=next_point_to_probe,
															BEST_SCORE=BEST_SCORE, BEST_PARAMS=BEST_PARAMS)

	# create optimised model and save
	model.params = BEST_PARAMS
	model.train()

	outname = id + '_model.pkl'
	pickle.dump(model, open(outname, "wb"))
	generic.save_var_models(model, outname, args['cv_fold'])
	print('Optimised model saved in ', outname)

	return dataset, BEST_SCORE

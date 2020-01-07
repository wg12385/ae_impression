import bayes_opt as bayes
from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction
from bayes_opt.observer import JSONLogger
from bayes_opt.event import Events
from bayes_opt.util import load_logs

import pickle
import numpy as np

import copy

from ml.hyperparameter_tuning import HPS_generic as generic

def load_logs(args):
	to_load = []
	try:
		with open(args['logfile'], 'r') as f:
			switch = False
			for line in f:
				if 'START' in line:
					switch = True
					print('SWITCH')
				elif switch:
					items = line.split()
					if items[0] == 'i':
						headers = items[1:]

						headers = [x.strip(' ') for x in headers]
						print(headers)

					else:
						params = {}
						print(items)
						for i in range(2, len(headers)):
							print(headers[i-1])
							if headers[i-1] == 'Mins':
								continue
							params[headers[i-1]] = items[i]
						print(params)
						score = items[1]

					to_load.append([params, score])
	except exception as e:
		print(e)
		return []

	return to_load


def gaussian_search(dataset, args):

	# determine whether log dictionary was provided
	if len(args['param_logs']) == 0:
		check_logs = False
	else:
		check_logs = True

	pbounds = {}
	for param in args['param_ranges'].keys():
		if args['param_logs'][param] == 'no':
			continue
		#pbounds[param] = (args['param_ranges'][param][0], args['param_ranges'][param][1])
		pbounds[param] = (0, 1)

	optimizer = BayesianOptimization(
		f=None,
		pbounds=pbounds,
		verbose=0, # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
		random_state=None
	)
	utility = UtilityFunction(kind="ucb", kappa=args['kappa'], xi=args['xi'])

	if args['load']:
		print('Loading previous data. . .')
		to_load = load_logs(args)
		for log in to_load:
			optimiser.register(params=log[0], target=(99999.99-logs[1])/99999.99)
		print('Loaded ', len(to_load), ' datapoints')

	generic.setup_logfile(args)

	if args['load']:
		for params, score in to_load:
			with open(args['logfile'], 'a') as f:
				iter = -1
				string = '{i:<10d}\t{score:<10.5f}'.format(i=iter, score=score)
				for param in args['param_list']:
					if args['param_logs'][param] == 'no':
						continue
					string = string + '\t{param:<20.4g}'.format(param=params[param])
				print(string, file=f)

	BEST_SCORE = 999.999
	BEST_PARAMS = {}
	searched = []
	force_random = False
	for e in range(int(args['epochs'])):

		if force_random or (args['random'] > 0 and e%args['random'] == 0):
			print('RANDOM')
			next_point_to_probe = {}
			for param in args['param_ranges'].keys():
				next_point_to_probe[param] = np.random.uniform(0, 1)
				#next_point_to_probe[param] = np.random.uniform(args['param_ranges'][param][0], args['param_ranges'][param][1])
			force_random = False
		else:
			next_point_to_probe = optimizer.suggest(utility)


		for param in args['param_ranges'].keys():
			if args['param_logs'][param] == 'no':
				continue
			next_point_to_probe[param] = (next_point_to_probe[param]*(args['param_ranges'][param][1]-args['param_ranges'][param][0])) + args['param_ranges'][param][0]

		if check_logs:
			for param in args['param_ranges'].keys():
				if args['param_logs'][param] == 'no':
					continue
				if 'log' in args['param_logs'][param]:
					next_point_to_probe[param] = 10**next_point_to_probe[param]

		score, BEST_SCORE, BEST_PARAMS = generic.HPS_iteration(e, dataset, args, next_point_to_probe=copy.copy(next_point_to_probe),
															BEST_SCORE=BEST_SCORE, BEST_PARAMS=BEST_PARAMS)

		searched.append(next_point_to_probe)
		try:
			optimizer.register(params=next_point_to_probe, target=(99999.99-score)/99999.99)
		except KeyError as e:
			print(e)
			print('Non-unique point generated, (increasing kappa value and) skipping (kappa = ',args['kappa'],')')
			args['kappa'] += 1.0
			utility = UtilityFunction(kind="ei", kappa=args['kappa'], xi=args['xi'])
			print('(forcing random point)')
			force_random = True


	outname = generic.save_models(dataset, BEST_PARAMS, args)
	print('Optimised model(s) saved in ', outname)

	return dataset, BEST_SCORE

























###

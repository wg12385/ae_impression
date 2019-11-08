from ml.models import FCHLmodel, KRRmodel, TFMmodel
import bayes_opt as bayes
from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction
from bayes_opt.observer import JSONLogger
from bayes_opt.event import Events
from bayes_opt.util import load_logs

import numpy as np


def gaussian_search(x, y, modelflag='KRR', featureflag='CMAT', id='test_model', logfile='gaussian.log',
 				param_ranges={}, param_logs={}, cv_steps=5, epochs=500, kappa=5, xi=0.1, random=0):

	# create model
	if modelflag == 'KRR':
		model = KRRmodel.KRRmodel(id, x, y, params={})
	elif modelflag == 'FCHL':
		model = FCHLmodel.FCHLmodel(id, x, y, params={})
	elif modelflag == 'TFM':
		model = TFMmodel.TFMmodel(id, x, y, params={})


	# determine whether log dictionary was provided
	if len(param_logs) == 0:
		check_logs = False
	else:
		check_logs = True

	pbounds = {}
	for param in param_ranges.keys():
		pbounds[param] = (param_ranges[param][0], param_ranges[param][1])

	optimizer = BayesianOptimization(
		f=None,
		pbounds=pbounds,
		verbose=0, # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
		random_state=None
	)
	utility = UtilityFunction(kind="ucb", kappa=kappa, xi=xi)

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



	BEST_SCORE = 999.999
	BEST_PARAMS = {}
	for _ in range(epochs):

		if random > 0 and _%random == 0:
			next_point_to_probe = {}
			for param in param_ranges.keys():
				next_point_to_probe[param] = np.random.uniform(param_ranges[param][0], param_ranges[param][1])

		else:
			next_point_to_probe = optimizer.suggest(utility)

		if check_logs:
			for param in param_ranges.keys():
				if param_logs[param] == True:
					next_point_to_probe[param] = 10**next_point_to_probe[param]

		model.params = next_point_to_probe

		y_pred = model.cv_predict(cv_steps)

		score = np.mean(np.absolute(y_pred - y))

		with open(logfile, 'a') as f:
			string = '{i:<10d}\t{score:<10.5f}'.format(i=p, score=score)
			for param in next_point_to_probe.keys():
				string = string + '\t{param:<15.4g}'.format(next_point_to_probe[param])

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

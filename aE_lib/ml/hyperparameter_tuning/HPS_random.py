
from ml.models import FCHLmodel, KRRmodel, TFMmodel
import numpy as np
import itertools
import pickle
import sys



def random_search(x, y, modelflag='KRR', featureflag='CMAT', targetflag='CCS', id='test_model', logfile='RANDOM.log',
 				param_ranges={}, param_logs={}, cv_steps=5, epochs=500):

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


	strings = []
	strings.append('HPS RANDOM SEARCH')
	strings.append(modelflag + '   ' + featureflag + '   ' + targetflag)
	for param in param_ranges.keys():
		strings.append('{param:<10s}: {low:>10.4g}  <--->  {high:<10.4g}  \t[{log:<s}]'.format(param=param,
																	low=param_ranges[param][0],
																	high=param_ranges[param][1]),
																	log=param_logs[param])
	strings.append('')
	strings.append('START')
	string = '{i:<10s}\t{score:<10s}'.format(i='i', score='SCORE')
	for param in param_ranges.keys():
		string = string + '\t{param:<15s}'.format(param=param)
	strings.append(string)
	with open(logfile, 'w') as f:
		for string in strings:
			print(string, file=f)

	BEST_SCORE = 999.999
	BEST_PARAMS = {}
	next_point_to_probe = {}

	for _ in range(epochs):
		for param in param_ranges.keys():
			next_point_to_probe[param] = np.random.uniform(param_ranges[param][0], param_ranges[param][1])

		if check_logs:
			for param in param_ranges.keys():
				if param_logs[param] == 'log':
					next_point_to_probe[param] = 10**next_point_to_probe[param]

		model.params = next_point_to_probe

		y_pred = model.cv_predict(cv_steps)

		score = np.mean(np.absolute(y_pred - y))

		with open(logfile, 'a') as f:
			string = '{i:<10d}\t{score:<10.5f}'.format(i=_, score=score)
			for param in next_point_to_probe.keys():
				string = string + '\t{param:<15.4g}'.format(param=next_point_to_probe[param])
			print(string, file=f)

		if score < BEST_SCORE:
			BEST_SCORE = score
			BEST_PARAMS = next_point_to_probe
			print('score = ', score)

	# create optimised model and save
	model.params = BEST_PARAMS
	model.train()

	outname = id + '_model.pkl'
	pickle.dump(model, open(outname, "wb"))

	print('Optimised model saved in ', outname)

	return score


def full_random_search(dataset, modelflag='KRR', featureflag='CMAT', targetflag='CCS', id='test_model', logfile='RANDOM.log',
 				param_ranges={}, param_logs={}, cv_steps=5, epochs=500):

	# determine whether log dictionary was provided
	if len(param_logs) == 0:
		check_logs = False
	else:
		check_logs = True


	strings = []
	strings.append('HPS RANDOM SEARCH')
	strings.append(modelflag + '   ' + featureflag + '   ' + targetflag)
	for param in param_ranges.keys():
		strings.append('{param:<10s}: {low:>10.4g}  <--->  {high:<10.4g}'.format(param=param,
																	low=param_ranges[param][0],
																	high=param_ranges[param][1]))
	strings.append('')
	strings.append('START')
	string = '{i:<10s}\t{score:<10s}'.format(i='i', score='SCORE')
	for param in param_ranges.keys():
		string = string + '\t{param:<15s}'.format(param=param)
	strings.append(string)
	with open(logfile, 'w') as f:
		for string in strings:
			print(string, file=f)

	BEST_SCORE = 999.999
	BEST_PARAMS = {}
	next_point_to_probe = {}

	for _ in range(epochs):
		for param in param_ranges.keys():
			next_point_to_probe[param] = np.random.uniform(param_ranges[param][0], param_ranges[param][1])

		if check_logs:
			for param in param_ranges.keys():
				if param_logs[param] == 'log':
					next_point_to_probe[param] = 10**next_point_to_probe[param]



		dataset.get_features_frommols(featureflag, targetflag, params=next_point_to_probe)

		assert len(dataset.x) > 0
		assert len(dataset.y) > 0

		# create model
		if modelflag == 'KRR':
			model = KRRmodel.KRRmodel(id, dataset.x, dataset.y, params=next_point_to_probe)
		elif modelflag == 'FCHL':
			model = FCHLmodel.FCHLmodel(id, dataset.x, dataset.y, params=next_point_to_probe)
		elif modelflag == 'TFM':
			model = TFMmodel.TFMmodel(id, dataset.x, dataset.y, params=next_point_to_probe)

		y_pred = model.cv_predict(cv_steps)

		score = np.mean(np.absolute(y_pred - dataset.y))

		with open(logfile, 'a') as f:
			string = '{i:<10d}\t{score:<10.5f}'.format(i=_, score=score)
			for param in next_point_to_probe.keys():
				string = string + '\t{param:<15.4g}'.format(param=next_point_to_probe[param])
			print(string, file=f)

		if score < BEST_SCORE:
			BEST_SCORE = score
			BEST_PARAMS = next_point_to_probe
			print('score = ', score)

	# create optimised model and save
	model.params = BEST_PARAMS
	model.train()

	outname = id + '_model.pkl'
	pickle.dump(model, open(outname, "wb"))

	print('Optimised model saved in ', outname)

	return dataset, BEST_SCORE

from ml.models import FCHLmodel, KRRmodel, NNmodel, TFMmodel
from sklearn.model_selection import KFold
import pickle
import numpy as np
import copy

def setup_logfile(args):
	strings = []
	strings.append('HPS ' + args['searchflag'] + ' SEARCH')
	strings.append(args['modelflag'] + '   ' + args['featureflag'] + '   ' + args['targetflag'])
	for param in args['param_ranges'].keys():
		strings.append('{param:<10s}: {low:>10.4g}  <--->  {high:<10.4g}'.format(param=param,
																	low=args['param_ranges'][param][0],
																	high=args['param_ranges'][param][1]))
	strings.append('')
	strings.append('START')
	string = '{i:<10s}\t{score:<10s}'.format(i='i', score='SCORE')
	for param in args['param_ranges'].keys():
		string = string + '\t{param:<20s}'.format(param=param)
	strings.append(string)

	with open(args['logfile'], 'w') as f:
		for string in strings:
			print(string, file=f)

def HPS_iteration(iter, dataset, args, next_point_to_probe={}, BEST_SCORE=999, BEST_PARAMS={}):

	print('HPS iter --- params', next_point_to_probe)

	print('HPS_iteration. . .')


	args['max_size'] = 200
	if args['featureflag'] == 'BCAI' and iter == 0:
		dataset.get_features_frommols(args, params=next_point_to_probe)
	elif args['feature_optimisation'] == 'True':
		dataset.get_features_frommols(args, params=next_point_to_probe)

	assert len(dataset.x) > 0
	assert len(dataset.y) > 0

	# create model
	if args['modelflag'] == 'KRR':
		model = KRRmodel.KRRmodel(id, dataset.x, dataset.y, params=next_point_to_probe, model_args=args)
	elif args['modelflag'] == 'FCHL':
		model = FCHLmodel.FCHLmodel(id, dataset.x, dataset.y, params=next_point_to_probe, model_args=args)
	elif args['modelflag'] == 'NN':
		model = NNmodel.NNmodel(id, dataset.x, dataset.y, params=next_point_to_probe, model_args=args)
	elif args['modelflag'] == 'TFM':
		model = TFMmodel.TFMmodel(id, dataset.x, dataset.y, params=next_point_to_probe, model_args=args)

	y_pred = model.cv_predict(args['cv_steps'])

	score = np.mean(np.absolute(y_pred - dataset.y))
	if score > 9999.99 or np.isnan(score):
		score = 999.99

	with open(args['logfile'], 'a') as f:
		string = '{i:<10d}\t{score:<10.5f}'.format(i=iter, score=score)
		for param in next_point_to_probe.keys():
			string = string + '\t{param:<20.4g}'.format(param=next_point_to_probe[param])
		print(string, file=f)

	if score < BEST_SCORE:
		BEST_SCORE = score
		BEST_PARAMS = next_point_to_probe
		print('BEST_SCORE = ', score)
	else:
		print('score = ', score)

	return score, BEST_SCORE, BEST_PARAMS

def save_models(dataset, BEST_PARAMS, args):

	if args['modelflag'] == 'KRR':
		model = KRRmodel.KRRmodel(id, dataset.x, dataset.y, params=BEST_PARAMS, model_args=args)
	elif args['modelflag'] == 'FCHL':
		model = FCHLmodel.FCHLmodel(id, dataset.x, dataset.y, params=BEST_PARAMS, model_args=args)
	elif args['modelflag'] == 'NN':
		model = NNmodel.NNmodel(id, dataset.x, dataset.y, params=BEST_PARAMS, model_args=args)
	elif args['modelflag'] == 'TFM':
		model = TFMmodel.TFMmodel(id, dataset.x, dataset.y, params=BEST_PARAMS, model_args=args)

	model.train()

	outname = args['output_dir'] +  args['modelflag'] + '_' + args['targetflag'] + '_' + args['featureflag'] + '_' + args['searchflag'] + '_model.pkl'

	pickle.dump(model, open(outname, "wb"))

	kf = KFold(n_splits=args['cv_steps'])
	kf.get_n_splits(model.train_x)

	i = 0
	for train_index, test_index in kf.split(model.train_x):
		i += 1

		tmp_model = copy.deepcopy(model)

		tmp_model.train_x = np.asarray(model.train_x)[train_index]
		tmp_model.train_y = np.asarray(model.train_y)[train_index]
		tmp_model.train()

		outfile = args['output_dir'] + outname.split('/')[-1].split('.')[0] + '_' + str(i) + '.pkl'
		pickle.dump(model, open(outfile, "wb"))

	return outname












#

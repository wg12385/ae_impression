from sklearn.model_selection import KFold
import pickle
import numpy as np
import copy
import time

def setup_logfile(args):
	strings = []
	strings.append('HPS ' + args['searchflag'] + ' SEARCH')
	strings.append(args['modelflag'] + '   ' + args['featureflag'] + '   ' + args['targetflag'])
	for param in args['param_list']:
		if args['param_logs'][param] == 'no':
			continue
		strings.append('{param:<10s}: {low:>10.4g}  <--->  {high:<10.4g}'.format(param=param,
																	low=args['param_ranges'][param][0],
																	high=args['param_ranges'][param][1]))
	strings.append('')
	strings.append('START')
	string = '{i:<10s}\t{score:<10s}'.format(i='i', score='SCORE')
	for param in args['param_list']:
		if args['param_logs'][param] == 'no':
			continue
		string = string + '\t{param:<20s}'.format(param=param)
	string = string + '\t{time:<10s}'.format(time='Mins')
	strings.append(string)

	if args['logfile'] == "":
		args['logfile'] = args['modelflag'] + '_' + args['featureflag'] + '_' + args['targetflag'] + args['searchflag'] + '.log'

	with open(args['logfile'], 'w') as f:
		for string in strings:
			print(string, file=f)

def HPS_iteration(iter, dataset, args, next_point_to_probe={}, BEST_SCORE=100000.00, BEST_PARAMS={}):

	print('HPS_iteration. . .')
	time0 = time.time()

	print('Params: ', next_point_to_probe)

	args['max_size'] = 200
	if args['featureflag'] == 'BCAI' and iter == 0:
		dataset.get_features_frommols(args, params=next_point_to_probe)
	elif args['feature_optimisation'] == 'True':
		dataset.get_features_frommols(args, params=next_point_to_probe)

	assert len(dataset.x) > 0
	assert len(dataset.y) > 0

	# create model
	# yes i know this is super bad "practice"
	# if you can think of a better way of allowing QML code to not be screwed up by TF or PyTorch, etc then let me know
	# otherwise, fight me . . .
	if args['modelflag'] == 'KRR':
		from ml.models import KRRmodel
		model = KRRmodel.KRRmodel(id, dataset.x, dataset.y, params=next_point_to_probe, model_args=args)
	elif args['modelflag'] == 'FCHL':
		from ml.models import FCHLmodel
		model = FCHLmodel.FCHLmodel(id, dataset.x, dataset.y, params=next_point_to_probe, model_args=args)
	elif args['modelflag'] == 'NN':
		from ml.models import NNmodel
		model = NNmodel.NNmodel(id, dataset.x, dataset.y, params=next_point_to_probe, model_args=args)
	elif args['modelflag'] == 'TFM':
		from ml.models import TFMmodel
		model = TFMmodel.TFMmodel(id, dataset.x, dataset.y, params=next_point_to_probe, model_args=args)

	y_pred = model.cv_predict(args['cv_steps'])

	score = np.mean(np.absolute(y_pred - dataset.y))
	if score > 99999.99 or np.isnan(score):
		score = 99999.99

	time1 = time.time()

	with open(args['logfile'], 'a') as f:
		string = '{i:<10d}\t{score:<10.5f}'.format(i=iter, score=score)
		for param in args['param_list']:
			if args['param_logs'][param] == 'no':
				continue
			string = string + '\t{param:<20.4g}'.format(param=next_point_to_probe[param])

		string = string + '\t{time:<10.4f}'.format(time=(time1-time0)/60)
		print(string, file=f)

	if score < BEST_SCORE:
		BEST_SCORE = score
		BEST_PARAMS = next_point_to_probe
		print('BEST_SCORE = ', score)
	else:
		print('score = ', score)

	if BEST_PARAMS == {}:
		BEST_PARAMS = next_point_to_probe

	return score, BEST_SCORE, BEST_PARAMS

def save_models(dataset, BEST_PARAMS, args):

	# create model
	if args['modelflag'] == 'KRR':
		from ml.models import KRRmodel
		model = KRRmodel.KRRmodel(id, np.asarray(dataset.x), np.asarray(dataset.y), params=BEST_PARAMS, model_args=args)
	elif args['modelflag'] == 'FCHL':
		from ml.models import FCHLmodel
		model = FCHLmodel.FCHLmodel(id, np.asarray(dataset.x), np.asarray(dataset.y), params=BEST_PARAMS, model_args=args)
	elif args['modelflag'] == 'NN':
		from ml.models import NNmodel
		model = NNmodel.NNmodel(id, dataset.x, dataset.y, params=BEST_PARAMS, model_args=args)
	elif args['modelflag'] == 'TFM':
		from ml.models import TFMmodel
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
		pickle.dump(tmp_model, open(outfile, "wb"))

	return outname












#

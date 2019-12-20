def HPS(dataset, args):

	if args['searchflag'] == 'grid':
		from .hyperparameter_tuning import HPS_grid
		dset, score = HPS_grid.grid_search(dataset, args)

	elif args['searchflag'] == 'gaussian':
		from .hyperparameter_tuning import HPS_gaussian
		dset, score = HPS_gaussian.gaussian_search(dataset, args)

	elif args['searchflag'] == 'random':
		from .hyperparameter_tuning import HPS_random
		dset, score = HPS_random.random_search(dataset, args)

	print('Model optimised, score = ', score)

	return dset, score

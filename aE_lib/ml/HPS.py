from .hyperparameter_tuning import HPS_gaussian, HPS_grid, HPS_random

def HPS(dataset, args):

	if args['searchflag'] == 'grid':
		dset, score = HPS_grid.grid_search(dataset, args)

	elif args['searchflag'] == 'gaussian':
		dset, score = HPS_gaussian.gaussian_search(dataset, args)

	elif args['searchflag'] == 'random':
		dset, score = HPS_random.random_search(dataset, args)

	print('Model optimised, score = ', score)

	return dset, score

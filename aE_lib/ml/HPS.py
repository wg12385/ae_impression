from .hyperparameter_tuning import HPS_gaussian, HPS_grid, HPS_random


def HPS(dataset, modelflag='FCHL', featureflag='FCHL', targetflag='CCS', searchflag='grid',
 				id='test_model', logfile='HPS.log',
				param_ranges={}, param_logs={}, grid_density=10, cv_steps=5,
				epochs=500, kappa=5, xi=0.1, random=0):

	if searchflag == 'grid':
		score = HPS_grid.grid_search(dataset.x, dataset.y, modelflag, featureflag, targetflag, id, logfile,
												param_ranges, param_logs, grid_density,
												int(cv_steps))

	elif searchflag == 'gaussian':
		score = HPS_gaussian.gaussian_search(dataset.x, dataset.y, modelflag, featureflag, targetflag, id, logfile,
												param_ranges, param_logs, int(cv_steps), int(epochs),
												kappa, xi, random)

	elif searchflag == 'random':
		score = HPS_random.random_search(dataset.x, dataset.y, modelflag, featureflag, targetflag, id, logfile,
												param_ranges, param_logs, int(cv_steps), int(epochs))

	print('Model optimised, score = ', score)


def Full_HPS(dataset, modelflag='FCHL', featureflag='FCHL', targetflag='CCS', searchflag='grid',
 				id='test_model', logfile='HPS.log',
				param_ranges={}, param_logs={}, grid_density=10, cv_steps=5,
				epochs=500, kappa=5, xi=0.1, random=0):


	if searchflag == 'grid':
		dset, score = HPS_grid.full_grid_search(dataset, modelflag, featureflag, targetflag, id, logfile,
												param_ranges, param_logs, grid_density, int(cv_steps))

	elif searchflag == 'gaussian':
		dset, score = HPS_gaussian.full_gaussian_search(dataset, modelflag, featureflag, targetflag, id, logfile,
												param_ranges, param_logs, int(cv_steps), int(epochs),
												kappa, xi, random)

	elif searchflag == 'random':
		dset, score = HPS_random.full_random_search(dataset, modelflag, featureflag, targetflag, id, logfile,
												param_ranges, param_logs, int(cv_steps), int(epochs))

	print('Model optimised, score = ', score)

	return dset, score

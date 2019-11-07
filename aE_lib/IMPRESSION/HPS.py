from .hyperparameter_tuning import HPS_gaussian, HPS_grid, HPS_random


def HPS(x, y, modelflag='FCHL', featureflag='FCHL', searchflag='grid',
 				id='test_model', logfile='HPS.log',
				param_ranges={}, param_logs={}, grid_density=10, cv_steps=5,
				epochs=500, kappa=5, xi=0.1, random=0):

	if searchflag == 'grid':
		score = HPS_grid.grid_search(x, y, modelflag, featureflag, id, logfile,
												param_ranges, param_logs, grid_density,
												cv_steps)

	elif searchflag == 'gaussian':
		score = HPS_gaussian.gaussian_search(x, y, modelflag, featureflag, id, logfile,
												param_ranges, param_logs, cv_steps, epochs,
												kappa, xi, random)

	elif searchflag == 'random':
		score = HPS_random.random_search(x, y, modelflag, featureflag, id, logfile,
												param_ranges, param_logs, cv_steps, epochs)

	print('Model optimised, score = ', score)

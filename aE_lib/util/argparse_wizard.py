
import glob
from util.flag_handler import hdl_targetflag, flag_combos, paramdict

def run_wizard(args):

	# code to run user through selection of input arguments for selected command
	if args['Command'] == 'setup_train':
		# Training set ##############################################################################
		check = False
		while not check:
			if args['training_set'].split('.')[-1] == 'pkl' or args['training_set'].split('.')[-1] == 'csv':
				try:
					a = open(args['training_set'], 'r')
					check = True
				except Exception as e:
					print(e)

			else:
				files = glob.glob(args['training_set'])
				try:
					a = open(files[0], 'r')
					check = True
				except Exception as e:
					print(e)

			if not check:
				args['training_set'] = input("Training set: \n")

		# Store datasets ##############################################################################
		check = False
		while not check:
			decision = input("Do you want to store dataset as pickle after creation ? [y]/n: \n")
			if len(decision) == 0:
				args['store_datasets'] = 'True'
				check = True
			elif decision[0] in ['n', 'N']:
				args['store_datasets'] = 'False'
				check = True
			elif decision[0] in ['y', 'Y']:
				args['store_datasets'] = 'True'
				check = True
			else:
				check = False

		combocheck = False
		while not combocheck:
			# Model ##############################################################################
			check = False
			while not check:
				model = input("What type of model do you want ? KRR, FCHL, NN, TFM : \n")
				if model in ['KRR', 'FCHL', 'NN', 'TFM']:
					args['modelflag'] = model
					check = True
			if model == 'FCHL':
				args['featureflag'] = 'FCHL'
				check = True
				combocheck = True
			else:
				# Feature ##############################################################################
				check = False
				while not check:
					feature = input("What type of input features do you want ? CMAT, aSLATM, FCHL, ACSF, BCAI : \n")
					if feature in ['CMAT', 'aSLATM', 'FCHL', 'ACSF', 'BCAI']:
						args['featureflag'] = feature
						check = True

				combocheck = flag_combos.check_combination(args['modelflag'], args['featureflag'])
				if not combocheck:
					print('Invalid combination of model and feature, try again . . .')

		# Target ##############################################################################
		check = False
		while not check:
			target = input("What is the target parameter ? XCS or nJXY : \n")
			param = hdl_targetflag.flag_to_target(target)
			if param == 0:
				print('Invalid parameter flag')
			else:
				args['targetflag'] = target
				check = True

		# Search method ##############################################################################
		check = False
		while not check:
			searchmethod = input("What search method should be used ? grid, gaussian, random : \n")
			if searchmethod in ['grid', 'gaussian', 'random']:
				args['searchmethod'] = searchmethod
				check = True

		# Feature optimisation ##############################################################################
		check = False
		while not check:
			feature_opt = input("Do you want to include feature parameters in optimisation ? : [y]/n \n")
			if len(feature_opt) == 0 or feature_opt[0] in ['Y', 'y']:
				args['feature_optimisation'] = 'True'
				check = True
			elif feature_opt[0] in ['N', 'n']:
				args['feature_optimisation'] = 'False'
				check = True

			# Feature file ##############################################################################
			if args['feature_optimisation'] == 'False':
				check = False
				while not check:
					file = input("File containing pre-made features dataset object: \n")
					try:
						a = open(file, 'r')
						check = True
					except Exception as e:
						print(e)

		args['param_ranges'], args['param_logs'] = paramdict.construct_param_dict(args['modelflag'], args['featureflag'], args['targetflag'])

		# Parameters ##############################################################################
		for param in args['param_ranges'].keys():
			check = False
			IP = input("Select range for parameter (min, max, log) {param:<10s}: default = {min:<10f}, {max:<10f}, {log:<10s} \n".format(param=param,
																												min=args['param_ranges'][param][0],
																												max=args['param_ranges'][param][1],
																												log=args['param_logs'][param]))
			if len(IP) == 0:
				check = True
			else:
				try:
					range = [float(IP.split(',')[0]), float(IP.split(',')[1])]
					log = IP.split(',')[2]

					args['param_ranges'][param] = range
					args['param_logs'][param] = log

					check = True

				except Exception as e:
					print(e)



		# grid density ##############################################################################
		check = False
		while not check:
			try:
				cv = input("Specify number of cross validation iterations: default = {0:<10f} \n".format(args['cv_steps']))
				if len(cv) == 0:
					check = True
				else:
					args['cv_steps'] = int(cv)
					check = True
			except Exception as e:
				print(e)


		if args['searchmethod'] == 'grid':
			# grid density ##############################################################################
			check = False
			while not check:
				try:
					grid = input("Specify grid density for parameters: default = {0:<10f} \n".format(args['grid_density']))
					if len(grid) == 0:
						check = True
					else:
						args['grid_density'] = int(grid)
						check = True
				except Exception as e:
					print(e)

		elif args['searchmethod'] == 'random':
			# epochs ##############################################################################
			check = False
			while not check:
				try:
					epochs = input("Specify number of epochs to run: default = {0:<10f} \n".format(args['epochs']))
					if len(epochs) == 0:
						check = True
					else:
						args['epochs'] = float(epochs)
						check = True
				except Exception as e:
					print(e)

		elif args['searchmethod'] == 'gaussian':
			# epochs ##############################################################################
			check = False
			while not check:
				try:
					epochs = input("Specify number of epochs to run: default = {0:<10f} \n".format(args['epochs']))
					if len(epochs) == 0:
						check = True
					else:
						args['epochs'] = float(epochs)
						check = True
				except Exception as e:
					print(e)

			# kappa ##############################################################################
			check = False
			while not check:
				try:
					kappa = input("Specify kappa value: default = {0:<10f} \n".format(args['kappa']))
					if len(kappa) == 0:
						check = True
					else:
						args['kappa'] = float(kappa)
						check = True
				except Exception as e:
					print(e)

			# xi ##############################################################################
			check = False
			while not check:
				try:
					xi = input("Specify xi value: default = {0:<10f} \n".format(args['xi']))
					if len(xi) == 0:
						check = True
					else:
						args['xi'] = float(xi)
						check = True
				except Exception as e:
					print(e)
			# random ##############################################################################
			check = False
			while not check:
				try:
					random = input("Specify frequency of random samples: default = {0:<10d} \n".format(args['random']))
					if len(random) == 0:
						check = True
					else:
						args['random'] = int(random)
						check = True
				except Exception as e:
					print(e)


	return args






















###

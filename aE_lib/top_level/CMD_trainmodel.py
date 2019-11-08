
from file_creation.HPC_submission import make_HPC_header
from ml.HPS import HPS
from molecule.dataset import dataset
from molecule.nmrmol import nmrmol
import sys

# creates a submission file that will run a hyper parameter search
def setup_trainmodel(args):

	# NEED TO MAKE GPU VERSION
	python_env = args.python_env
	system = args.system
	nodes = 1
	ppn = args.processors
	walltime = args.walltime
	mem = args.memory

	jobname = 'IMP_' + args.modelflag + '_' + args.featureflag + '_' + args.targetflag + '_' + args.searchmethod + '_HPS'
	submission_file = jobname + '.submit'
	header = make_HPC_header(jobname, system=args.system, nodes=1, ppn=args.processors, walltime=args.walltime, mem=args.memory)
	strings = []
	strings.append('')
	strings.append('shopt -s expand_aliases')
	strings.append('source ~/.bashrc')
	strings.append("source activate {env:<10s}".format(env=python_env))
	strings.append("IMPRESSION train --modelflag '{model:<s}'  --featureflag '{feature:<s}' \\".format(model=args.modelflag,
																									feature=args.featureflag))
	strings.append("				 --targetflag '{param:<s}'  --searchmethod '{search:<s}'  \\".format(param=args.targetflag,
																									search=args.searchmethod))
	strings.append("				 --param_ranges '{param:<s}'  \\".format(param=str(args.param_ranges)))
	strings.append("				 --param_logs '{param:<s}'  \\".format(param=str(args.param_logs)))
	strings.append("				 --cv_steps {cv_steps:<d} \\".format(cv_steps=args.cv_steps))

	if args.searchmethod == "grid":
		strings.append("				 --grid_density {grid_density:<d} \\".format(grid_density=args.grid_density))

	elif args.searchmethod == "gaussian":
		strings.append("				 --epochs {epochs:<d} \\".format(epochs=args.epochs))
		strings.append("				 --kappa {kappa:<f} \\".format(kappa=args.kappa))
		strings.append("				 --xi {xi:<f} \\".format(xi=args.xi))
		strings.append("				 --random {random:<d} \\".format(random=args.random))

	elif args.searchmethod == "random":
		strings.append("				 --epochs {epochs:<d} \\".format(epochs=args.epochs))


	with open(submission_file, 'w') as f:
		for string in header:
			print(string, file=f)
		for string in strings:
			print(string, file=f)


def trainmodel(args):

	print(args.training_set)
	sys.exit(0)


	if len(args.training_set) == 1:
		if args.training_set.split('.')[-1] == 'pkl':
			dset = pickle.load(open(args.training_set, 'rb'))
		elif args.training_set.split('.')[-1] == 'csv':
			dset = load_dataset_from_csv(args.training_set)
	else:
		files = glob.glob(args.training_set)
		dset = dataset()
		if args.store_datasets:
			dataset.get_mols(files, type='nmredata')
			dataset.get_features_frommols(files, args.featureflag, args.targetflag, args.cutoff)
			pickle.dump(dataset, open('training_set.pkl', 'wb'))
		else:
			dataset.get_features_fromfiles(files, args.featureflag, args.targetflag, args.cutoff)


	HPS(dataset.x, dataset.y, modelflag=args.modelflag, featureflag=args.featureflag, searchflag=args.searchmethod,
 				id='test_model', logfile=args.logfile, param_ranges=args.param_ranges,
				param_logs=args.param_logs, grid_density=args.grid_density, cv_steps=args.cv_steps,
				epochs=args.epochs, kappa=args.kappa, xi=args.xi, random=args.random)





###
